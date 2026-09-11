# HydroTuring entrypoint for wflow_sbm (Wflow.jl v1.0.4).
#
#     julia --project=/opt/wflow_sbm ht_adapter.jl --request /io/request.json
#
# The adapter is the package in src/WflowSbmAdapter.jl, so that its code, and the Wflow
# code it calls for a daily and an hourly case, is compiled when the image is built and
# not every time a container starts.
using WflowSbmAdapter

exit(
    try
        WflowSbmAdapter.main(ARGS)
    catch err
        showerror(stderr, err, catch_backtrace())
        println(stderr)
        1
    end,
)
