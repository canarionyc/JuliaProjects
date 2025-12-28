julia --version

juliaup status

julia +1.10 -e "println(Sys.BINDIR)"

julia +1.10 -e "import Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()"

julia +1.10 -e 'import Pkg; Pkg.activate("."); Pkg.update("ApproxFun")'