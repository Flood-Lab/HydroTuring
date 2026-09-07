"""SAC-SMA, Snow-17 and a gamma unit hydrograph, in the standard library.

A line-by-line port of the NWS Fortran in Upstream-Tech/SACSMA-SNOW17
(sacsma_source/sac/sac1.f, ex_sac1.f, duamel.f and snow19/*.f, commit
8737b92), kept in the Fortran's own variable names so it can be read
against the original. The port was checked against the f2py build of that
Fortran on the benchmark's synthetic records before it was trusted.

What is deliberately left out: the frozen-ground option (IFRZE=0 in the
wrapper), the observation-updating routines (UPDT19, ADJC19, AECO19: they
need observed snow), the snow depth and density diagnostics (SNDEPTH,
SNOWPACK, SNOWT, SNEW: they move no water), and the user-specified melt
factor curve (LMFV=0). Nothing that moves water is omitted.

Units: depths in mm, temperature in degC, time in the step the caller
passes: SAC-SMA takes the step in days, Snow-17 in whole hours as the
original does.
"""

from __future__ import annotations

import math

# ----------------------------------------------------------------------------
# SAC-SMA (sac1.f), one time interval
# ----------------------------------------------------------------------------

SAC_PARAMS = ("uztwm", "uzfwm", "uzk", "pctim", "adimp", "riva", "zperc", "rexp",
              "lztwm", "lzfsm", "lzfpm", "lzsk", "lzpk", "pfree", "side", "rserv")


class SacState:
    __slots__ = ("uztwc", "uzfwc", "lztwc", "lzfsc", "lzfpc", "adimc")

    def __init__(self, uztwc=0.0, uzfwc=0.0, lztwc=0.0, lzfsc=0.0, lzfpc=0.0, adimc=0.0):
        self.uztwc, self.uzfwc, self.lztwc = uztwc, uzfwc, lztwc
        self.lzfsc, self.lzfpc, self.adimc = lzfsc, lzfpc, adimc


def sac1(dt: float, pxv: float, ep: float, p: dict, s: SacState) -> dict:
    """One SAC-SMA interval of `dt` days with `pxv` mm of rain plus melt and
    `ep` mm of demand. Mutates the state; returns the fluxes of sac1.f:
    tci (channel inflow), roimp, sdro, ssur, sif, bfs, bfp, tet, plus bfncc,
    the baseflow lost to deep groundwater when side > 0, which the Fortran
    computes and drops."""
    uztwm, uzfwm, uzk, pctim, adimp, riva = p["uztwm"], p["uzfwm"], p["uzk"], p["pctim"], p["adimp"], p["riva"]
    zperc, rexp, lztwm, lzfsm, lzfpm = p["zperc"], p["rexp"], p["lztwm"], p["lzfsm"], p["lzfpm"]
    lzsk, lzpk, pfree, side, rserv = p["lzsk"], p["lzpk"], p["pfree"], p["side"], p["rserv"]
    uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc = s.uztwc, s.uzfwc, s.lztwc, s.lzfsc, s.lzfpc, s.adimc

    # Evapotranspiration from the upper zone.
    edmnd = ep
    e1 = edmnd * (uztwc / uztwm)
    red = edmnd - e1
    uztwc -= e1
    e2 = 0.0
    if uztwc < 0.0:
        e1 += uztwc
        uztwc = 0.0
        red = edmnd - e1
        if uzfwc >= red:
            e2 = red
            uzfwc -= e2
            red = 0.0
        else:
            e2 = uzfwc
            uzfwc = 0.0
            red -= e2
    else:
        if (uztwc / uztwm) < (uzfwc / uzfwm):
            uzrat = (uztwc + uzfwc) / (uztwm + uzfwm)
            uztwc = uztwm * uzrat
            uzfwc = uzfwm * uzrat
    if uztwc < 0.00001:
        uztwc = 0.0
    if uzfwc < 0.00001:
        uzfwc = 0.0

    # Evapotranspiration from the lower zone.
    e3 = red * (lztwc / (uztwm + lztwm))
    lztwc -= e3
    if lztwc < 0.0:
        e3 += lztwc
        lztwc = 0.0
    ratlzt = lztwc / lztwm
    saved = rserv * (lzfpm + lzfsm)
    ratlz = (lztwc + lzfpc + lzfsc - saved) / (lztwm + lzfpm + lzfsm - saved)
    if ratlzt < ratlz:
        delta = (ratlz - ratlzt) * lztwm
        lztwc += delta
        lzfsc -= delta
        if lzfsc < 0.0:
            lzfpc += lzfsc
            lzfsc = 0.0
    if lztwc < 0.00001:
        lztwc = 0.0

    # Evapotranspiration from the ADIMP area.
    e5 = e1 + (red + e2) * ((adimc - e1 - uztwc) / (uztwm + lztwm))
    adimc -= e5
    if adimc < 0.0:
        e5 += adimc
        adimc = 0.0
    e5 *= adimp

    # Percolation and runoff.
    twx = pxv + uztwc - uztwm
    if twx < 0.0:
        uztwc += pxv
        twx = 0.0
    else:
        uztwc = uztwm
    adimc += pxv - twx
    roimp = pxv * pctim

    sbf = ssur = sif = sperc = sdro = spbf = 0.0
    ninc = int(1.0 + 0.2 * (uzfwc + twx))
    dinc = (1.0 / ninc) * dt
    pinc = twx / ninc
    duz = 1.0 - (1.0 - uzk) ** dinc
    dlzp = 1.0 - (1.0 - lzpk) ** dinc
    dlzs = 1.0 - (1.0 - lzsk) ** dinc
    parea = 1.0 - adimp - pctim

    for _ in range(ninc):
        adsur = 0.0
        ratio = (adimc - uztwc) / lztwm
        if ratio < 0.0:
            ratio = 0.0
        addro = pinc * ratio ** 2

        bf = lzfpc * dlzp
        lzfpc -= bf
        if lzfpc <= 0.0001:
            bf += lzfpc
            lzfpc = 0.0
        sbf += bf
        spbf += bf
        bf = lzfsc * dlzs
        lzfsc -= bf
        if lzfsc <= 0.0001:
            bf += lzfsc
            lzfsc = 0.0
        sbf += bf

        if (pinc + uzfwc) > 0.01:
            percm = lzfpm * dlzp + lzfsm * dlzs
            perc = percm * (uzfwc / uzfwm)
            defr = 1.0 - (lztwc + lzfpc + lzfsc) / (lztwm + lzfpm + lzfsm)
            perc = perc * (1.0 + zperc * defr ** rexp)
            if perc >= uzfwc:
                perc = uzfwc
            uzfwc -= perc
            check = lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm
            if check > 0.0:
                perc -= check
                uzfwc += check
            sperc += perc

            delta = uzfwc * duz
            sif += delta
            uzfwc -= delta

            perct = perc * (1.0 - pfree)
            if (perct + lztwc) > lztwm:
                percf = perct + lztwc - lztwm
                lztwc = lztwm
            else:
                lztwc += perct
                percf = 0.0
            percf += perc * pfree
            if percf != 0.0:
                hpl = lzfpm / (lzfpm + lzfsm)
                ratlp = lzfpc / lzfpm
                ratls = lzfsc / lzfsm
                fracp = (hpl * 2.0 * (1.0 - ratlp)) / ((1.0 - ratlp) + (1.0 - ratls))
                if fracp > 1.0:
                    fracp = 1.0
                percp = percf * fracp
                percs = percf - percp
                lzfsc += percs
                if lzfsc > lzfsm:
                    percs = percs - lzfsc + lzfsm
                    lzfsc = lzfsm
                lzfpc += percf - percs
                if lzfpc > lzfpm:
                    excess = lzfpc - lzfpm
                    lztwc += excess
                    lzfpc = lzfpm

            if pinc != 0.0:
                if (pinc + uzfwc) > uzfwm:
                    sur = pinc + uzfwc - uzfwm
                    uzfwc = uzfwm
                    ssur += sur * parea
                    adsur = sur * (1.0 - addro / pinc)
                    ssur += adsur * adimp
                else:
                    uzfwc += pinc
        else:
            uzfwc += pinc

        adimc += pinc - addro - adsur
        if adimc > (uztwm + lztwm):
            addro += adimc - (uztwm + lztwm)
            adimc = uztwm + lztwm
        sdro += addro * adimp
        if adimc < 0.00001:
            adimc = 0.0

    eused = e1 + e2 + e3
    sif *= parea
    tbf = sbf * parea
    bfcc = tbf * (1.0 / (1.0 + side))
    bfp = spbf * parea / (1.0 + side)
    bfs = bfcc - bfp
    if bfs < 0.0:
        bfs = 0.0
    bfncc = tbf - bfcc
    tci = roimp + sdro + ssur + sif + bfcc
    e4 = (edmnd - eused) * riva
    tci -= e4
    if tci < 0.0:
        e4 += tci
        tci = 0.0
    eused *= parea
    tet = eused + e5 + e4
    if adimc < uztwc:
        adimc = uztwc

    s.uztwc, s.uzfwc, s.lztwc, s.lzfsc, s.lzfpc, s.adimc = uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc
    return {"tci": tci, "roimp": roimp, "sdro": sdro, "ssur": ssur, "sif": sif,
            "bfs": bfs, "bfp": bfp, "tet": tet, "bfncc": bfncc, "parea": parea}


def sac_storage(p: dict, s: SacState) -> tuple[float, float]:
    """Catchment-average water in the soil column and in the lower-zone free
    water, weighting each store by the area it describes (the pervious area
    for the five zone stores, ADIMP for the additional impervious store)."""
    parea = 1.0 - p["adimp"] - p["pctim"]
    soil = parea * (s.uztwc + s.uzfwc + s.lztwc) + p["adimp"] * s.adimc
    lower = parea * (s.lzfsc + s.lzfpc)
    return soil, lower


# ----------------------------------------------------------------------------
# Snow-17 (exsnow19.f, PACK19.f, MELT19.f, AESC19.f, ROUT19.f, ZERO19.f)
# ----------------------------------------------------------------------------

SNOW_PARAMS = ("scf", "mfmax", "mfmin", "uadj", "si", "nmf", "tipm", "mbase", "pxtemp", "plwhc", "daygm")
JULDAY = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)


class SnowState:
    """The SNCO19 carryover: WE, NEGHS, LIQW, TINDEX, ACCMAX, SB, SBAESC, SBWS,
    STORGE, AEADJ, EXLAG(NEXLAG), plus TPREV."""

    def __init__(self, nexlag: int):
        self.we = self.neghs = self.liqw = self.tindex = self.accmax = 0.0
        self.sb = self.sbaesc = self.sbws = self.storge = self.aeadj = 0.0
        self.exlag = [0.0] * nexlag
        self.tprev = 0.0

    def zero(self):
        self.we = self.neghs = self.liqw = self.tindex = self.accmax = 0.0
        self.sb = self.sbaesc = self.sbws = self.storge = self.aeadj = 0.0
        for i in range(len(self.exlag)):
            self.exlag[i] = 0.0

    def total(self) -> float:
        return self.we + self.liqw + sum(self.exlag) + self.storge


def _day_from_march21(year: int, month: int, day: int) -> int:
    kda = JULDAY[month - 1] + day
    nda = 365
    if year % 4 == 0 and month >= 3:
        kda += 1
        nda += 1
    i0 = JULDAY[2] + 21
    i1 = JULDAY[month - 1] + day
    return i1 - i0 if kda >= i0 else nda - (i0 - i1)


def melt19(idn: int, alat: float, ta: float, mfmax: float, mfmin: float, mbase: float,
           tindex: float, tipm: float, nmf: float) -> tuple[float, float, float]:
    """Surface melt at full cover under non-rain conditions, the negative heat
    exchange, and the updated temperature index (MELT19, LMFV=0)."""
    diff = mfmax - mfmin
    dayn = float(idn)
    if alat < 54.0:
        mf = math.sin(dayn * 2.0 * 3.1416 / 366.0) * diff * 0.5 + (mfmax + mfmin) * 0.5
    else:
        if idn >= 275:
            x = (dayn - 275.0) / (458.0 - 275.0)
        elif idn >= 92:
            x = (275.0 - dayn) / (275.0 - 92.0)
        else:
            x = (91.0 + dayn) / 183.0
        xx = math.sin(dayn * 2.0 * 3.1416 / 366.0) * 0.5 + 0.5
        if x <= 0.48:
            adjmf = 0.0
        elif x >= 0.70:
            adjmf = 1.0
        else:
            adjmf = (x - 0.48) / (0.70 - 0.48)
        mf = (xx * adjmf) * diff + mfmin
    ratio = mf / mfmax
    tmx = ta - mbase
    if tmx < 0.0:
        tmx = 0.0
    tsur = ta if ta <= 0.0 else 0.0
    tnmx = tindex - tsur
    cnhs = ratio * nmf * tnmx
    tindex = tindex + tipm * (ta - tindex)
    if tindex > 0.0:
        tindex = 0.0
    melt = mf * tmx if tmx > 0.0 else 0.0
    return melt, cnhs, tindex


def aesc19(st: SnowState, si: float, adc: list[float], snof: float) -> float:
    twe = st.we + st.liqw
    if twe > st.accmax:
        st.accmax = twe
    if twe >= st.aeadj:
        st.aeadj = 0.0
    ai = st.accmax
    if st.accmax > si:
        ai = si
    if st.aeadj > 0.0:
        ai = st.aeadj
    if twe >= ai:
        st.sb = twe
        st.sbws = twe
        aesc = 1.0
    elif twe <= st.sb + 1e-6:  # the Fortran adds an unset `tiny` fudge; a small epsilon here
        r = (twe / ai) * 10.0 + 1.0
        n = int(r)
        r -= n
        aesc = adc[n - 1] + (adc[n] - adc[n - 1]) * r
        if aesc > 1.0:
            aesc = 1.0
        st.sb = twe + snof
        st.sbws = twe
        st.sbaesc = aesc
    elif twe >= st.sbws:
        aesc = 1.0
    else:
        aesc = st.sbaesc + (1.0 - st.sbaesc) * ((twe - st.sb) / (st.sbws - st.sb))
    if aesc < 0.05:
        aesc = 0.05
    if aesc > 1.0:
        aesc = 1.0
    return aesc


def rout19(it: int, excess: float, we: float, aesc: float, st: SnowState) -> float:
    """Lag and attenuate excess water through the pack (ROUT19). Returns the
    pack outflow for the interval; mutates storge and exlag."""
    fit = float(it)
    packro = 0.0
    cl = 0.03 * fit / 6.0
    exlag = st.exlag
    if excess != 0.0:
        if excess < 0.1 or we < 1.0:
            exlag[0] += excess
        else:
            n = int((excess * 4.0) ** 0.3 + 0.5)
            if n == 0:
                n = 1
            fn = float(n)
            for i in range(1, n + 1):
                fi = float(i)
                term = cl * we * fn / (excess * (fi - 0.5))
                if term > 150.0:
                    term = 150.0
                flag = 5.33 * (1.0 - math.exp(-term))
                l2 = int((flag + fit) / fit + 1.0)
                l1 = l2 - 1
                endl1 = l1 * it
                por2 = (flag + fit - endl1) / fit
                por1 = 1.0 - por2
                exlag[l2 - 1] += por2 * excess / fn
                exlag[l1 - 1] += por1 * excess / fn
    if (st.storge + exlag[0]) != 0.0:
        if (st.storge + exlag[0]) < 0.1:
            packro = st.storge + exlag[0]
            st.storge = 0.0
        else:
            el = exlag[0] / fit
            els = el / (25.4 * aesc)
            wes = we / (25.4 * aesc)
            term = 500.0 * els / (wes ** 1.3)
            if term > 150.0:
                term = 150.0
            r1 = 1.0 / (5.0 * math.exp(-term) + 1.0)
            for _ in range(it):
                os_ = (st.storge + el) * r1
                packro += os_
                st.storge = st.storge + el - os_
            if st.storge <= 0.001:
                packro += st.storge
                st.storge = 0.0
    for i in range(1, len(exlag)):
        exlag[i - 1] = exlag[i]
    exlag[-1] = 0.0
    return packro


def snow17(idt: int, year: int, month: int, day: int, pcp: float, tmp: float,
           alat: float, p: dict, pa: float, adc: list[float],
           st: SnowState) -> tuple[float, float]:
    """One Snow-17 interval of `idt` whole hours (EXSNOW19 + PACK19 with
    NDT = 1, no updating, no user melt curve). Returns (rain plus melt
    reaching the soil, snowfall). Mutates the state. Written in the
    Fortran's branch order; the labels in comments are PACK19's."""
    idn = _day_from_march21(year, month, day)
    mfmax = p["mfmax"] * idt / 6.0
    mfmin = p["mfmin"] * idt / 6.0
    nmf = p["nmf"] * idt / 6.0
    uadj = p["uadj"] * idt / 6.0
    pgm = p["daygm"] * idt / 24.0
    tipm = 1.0 - (1.0 - p["tipm"]) ** (idt / 6.0)
    si, scf, mbase, pxtemp, plwhc = p["si"], p["scf"], p["mbase"], p["pxtemp"], p["plwhc"]
    snof = 0.0  # SNUP19's SNOF is never set by the wrapper: zero

    ta = tmp
    itpx = idt
    fitpx = float(itpx)
    gm = pgm
    sfnew = 1.5 * fitpx
    rfmin = 0.25 * fitpx
    sbci = 0.0612 * fitpx
    pxi = pcp

    if pxi == 0.0 and st.we == 0.0:
        # 160: no snow cover, no new snowfall.
        return _end_of_period(st, ta, pxi, 0.0, si, adc, snof)

    sfall = cnhspx = rain = rainm = 0.0
    if pxi != 0.0:
        tpx = ta
        if tpx > pxtemp:
            fracs, fracr = 0.0, 1.0
        else:
            fracs, fracr = 1.0, 0.0
        if fracs != 0.0:
            ts = tpx if tpx <= 0.0 else 0.0
            sfall = pxi * fracs * scf
            if (st.we + st.liqw) >= st.sbws:
                st.sbws += 0.75 * sfall
                if sfall >= snof and st.sb > st.we + st.liqw:
                    st.sb = st.we + st.liqw
            elif sfall >= snof:
                st.sbws = st.we + st.liqw + 0.75 * sfall
            st.we += sfall
            if st.we + st.liqw >= 3.0 * st.sb:
                st.accmax = st.we + st.liqw
                st.aeadj = 0.0
            cnhspx = -ts * sfall / 160.0
            if sfall > sfnew:
                st.tindex = ts
        rain = pxi * fracr
        if st.we == 0.0:
            # 160: rain on bare ground.
            return _end_of_period(st, ta, pxi, sfall, si, adc, snof)
        tr = tpx if tpx >= 0.0 else 0.0
        rainm = 0.0125 * rain * tr

    # 110: melt at the ground-snow interface.
    if st.we <= gm:
        gmro = st.we + st.liqw
        melt = 0.0
        robg = rain
        rain = 0.0
        # 150: snow gone.
        packro = gmro + melt + sum(st.exlag) + st.storge + rain
        st.zero()
        return _end_of_period(st, ta, packro + robg, sfall, si, adc, snof)

    gmwlos = (gm / st.we) * st.liqw
    gmslos = gm
    pmelt, pcnhs, st.tindex = melt19(idn, alat, ta, mfmax, mfmin, mbase, st.tindex, tipm, nmf)
    cnhs = pcnhs
    if rain > rfmin:
        # 120: rain-on-snow energy balance.
        ea = 0.90 * 2.7489e8 * math.exp(-4278.63 / (ta + 242.792))
        tak = (ta + 273) * 0.01
        qn = sbci * (tak ** 4 - 55.55)
        qe = 8.5 * (ea - 6.11) * uadj
        qh = 7.5 * 0.000646 * pa * uadj * ta
        melt = qn + qe + qh + rainm
        if melt < 0.0:
            melt = 0.0
    else:
        melt = pmelt + rainm

    # 130: areal extent at the start of the interval.
    aesc = aesc19(st, si, adc, snof)
    if aesc != 1.0:
        melt *= aesc
        cnhs *= aesc
        gmwlos *= aesc
        gmslos *= aesc
        robg = (1.0 - aesc) * rain
        rain -= robg
    else:
        robg = 0.0
    if (cnhs + st.neghs) < 0.0:
        cnhs = -1.0 * st.neghs

    st.we -= gmslos
    st.liqw -= gmwlos
    gmro = gmslos + gmwlos

    if melt > 0.0 and melt >= st.we:
        # 150: the surface melt takes the whole pack.
        melt = st.we + st.liqw
        packro = gmro + melt + sum(st.exlag) + st.storge + rain
        st.zero()
        return _end_of_period(st, ta, packro + robg, sfall, si, adc, snof)
    if melt > 0.0:
        st.we -= melt

    # 137: heat and water balance of the cover.
    water = melt + rain
    heat = cnhs + cnhspx
    liqwmx = plwhc * st.we
    st.neghs += heat
    if st.neghs < 0.0:
        st.neghs = 0.0
    if st.neghs > 0.33 * st.we:
        st.neghs = 0.33 * st.we
    if (water + st.liqw) >= (liqwmx + st.neghs + plwhc * st.neghs):
        excess = water + st.liqw - liqwmx - st.neghs - plwhc * st.neghs
        st.liqw = liqwmx + plwhc * st.neghs
        st.we += st.neghs
        st.neghs = 0.0
    elif water >= st.neghs:
        st.liqw += water - st.neghs
        st.we += st.neghs
        st.neghs = 0.0
        excess = 0.0
    else:
        st.we += water
        st.neghs -= water
        excess = 0.0
    if st.neghs == 0.0:
        st.tindex = 0.0
    packro = rout19(itpx, excess, st.we, aesc, st) + gmro
    return _end_of_period(st, ta, packro + robg, sfall, si, adc, snof)


def _end_of_period(st: SnowState, ta: float, rm: float, sfall: float,
                   si: float, adc: list[float], snof: float) -> tuple[float, float]:
    """PACK19's close: the areal extent is recomputed from the conditions at
    the end of the period (which moves SB, SBWS and SBAESC), and TPREV is
    carried."""
    if st.total() != 0.0:
        aesc19(st, si, adc, snof)
    st.tprev = ta
    return rm, sfall


# ----------------------------------------------------------------------------
# Gamma unit hydrograph (duamel.f), with a real step
# ----------------------------------------------------------------------------

def gamma_uh(shape: float, scale_days: float, dt_days: float, max_len: int = 1000) -> list[float]:
    """DUAMEL's kernel: u(i) ∝ (i·dt/scale)^(shape−1)·exp(−i·dt/scale)/(Γ(shape)·scale),
    truncated where the log-ordinate drops below −8, normalised to one. The
    Fortran takes an integer step of days; this takes a real one so the
    hydrograph keeps its shape in days at any step."""
    if shape < 0.0:
        return [1.0]
    toc = math.log(math.gamma(shape) * scale_days)
    u = []
    for i in range(1, max_len + 1):
        top = i * dt_days / scale_days
        tor = (shape - 1.0) * math.log(top) - top - toc
        if tor > -8.0:
            u.append(math.exp(tor))
        elif i > 1:
            break
        else:
            u.append(0.0)
    total = sum(u) or 1e-5
    return [x / total for x in u]
