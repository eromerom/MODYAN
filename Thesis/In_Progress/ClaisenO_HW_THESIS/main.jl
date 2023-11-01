#!/usr/bin/env julia

include("gromacs_reader.jl")
include("sol_module.jl")
include("unk_module.jl")
include("misc_module.jl")
include("main_sol_module.jl")
include("main_unk_module.jl")
include("main_misc_module.jl")

println("Please enter 'SOL', 'UNK' or 'MISC' to start main program:")
user_input = readline()

if uppercase(user_input) == "SOL"
	main_sol()
elseif uppercase(user_input) == "UNK"
	main_unk()
elseif uppercase(user_input) == "MISC"
	main_misc()
else
	println("INVALID INPUT! Please enter either 'SOL', 'UNK' or 'MISC'.")
end

