
function show_help_sol()
	println("""
	┌─────────────────────────────────────────────────────────────────────────────────────┐
	│                             List of Available Functions                             │
	├─────────────────────────────────────────────────────────────────────────────────────┤
	│ 1: Count and update SOL molecule dictionaries.                                      │
	│ 2: Log SOL molecule coordinates.                                                    │
	│ 3: Log the center of mass for each SOL molecule.                                    │
	│ 4: Count and log the number of SOL molecules in a single cuboid.                    │
	│ 5: Count and log the number of SOL molecules in multiple cuboids.                   │
	│ 6: Calculate and log SOL molecule dipole moments in Cartesian coordinates.          │
	│ 7: Update and calculate the dipole moment library for SOL molecules in spherical    │
	│    coordinates.                                                                     │
	│ 8: Count and log the theta values for SOL molecules in multiple cuboids.            │
	│ 9: Generate RDF plots (SOL-SOL, single molecule).                                   │
	│ 10: Generate RDF plots (SOL-SOL, all molecules).                                    │
	│ help: Display the list of available functions.                                      │
	│ exit: Terminate the program.                                                        │
	└─────────────────────────────────────────────────────────────────────────────────────┘
	""")
end

function show_progress(stop_channel::Channel)
	print("Running")
	count = 0
	while true
		count += 1
		sleep(1)
		print(".")
		if count == 6  # Limit the number of dots printed
			count = 0
			print("\b" ^ 6)  # Erase the printed dots
		end

		if isready(stop_channel)
			break
		end
	end
	print("\nROUNTINE DONE! Type help or exit.\n")
end

function main_sol()
	println("╔════════════════════════════════════════════════════════════════════════════╗")
	println("║             WELCOME TO MODYAN (MOlecular DYnamics ANalyzer)                ║")
	println("╚════════════════════════════════════════════════════════════════════════════╝\n")
	println("This program is a Julia-based tool designed to analyze GROMACS files. The")
	println("format of the GROMACS file that this code can read is defined in the")
	println("subroutine called 'gromacs_reader.jl'. If you encounter errors, please")
	println("adjust the format reading loop in that subroutine.\n")
	println("To use this tool effectively, ensure that you have all the GROMACS file")
	println("snapshots in the current directory.\n")
	println("Copyright © April (2023) by Eduardo Romero-Montalvo and Alireza Sadeghifar,")
	println("University of British Columbia, Okanagan Campus.\n")
	println("This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.")
	println("This license allows you to share, copy, distribute, and transmit the work, as well as to remix,")
	println("transform, and build upon it for any purpose, provided the authors are given appropriate credit")
	println("and any derivative works are distributed under the same or similar license.\n")
	println("═════════════════════════════════════════════════════════════════════════════")
	println("Type 'help' for a list of available functions or 'exit' to quit the program.")
	println("═════════════════════════════════════════════════════════════════════════════\n")

	input_channel = Channel{String}(1)

	@async begin
		while true
			user_input = readline()
			put!(input_channel, user_input)
		end
	end

	while true
		input = take!(input_channel)


		if input == "exit"
			break
		elseif input == "help"
			show_help_sol()
		else

			function_number = tryparse(Int, input)

			if function_number == 1
				println("counting_sol() has been selected, no argument required.")

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				counting_sol()
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 2
				println("write_sol() has been selected, no argument required.")

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				write_sol()
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 3
				println("write_com_sol() has been selected, no arguments required.")

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				write_com_sol()
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 4
				println("write_sol_in_cuboid(iX, iY, iZ, fX, fY, fZ) has been selected.\n")
				println("Arguments are required.\n")
				println("Please enter them as iX   iY   iZ   fX   fY   fZ.\n")
				println("As a reference, the dimensions of the NVT simulation box are")
				print("0.0   0.0   0.0")

				gro_files = glob("*.gro")
				cmd1 = `tail -q -n 1 $gro_files`
				cmd2 = `uniq`

				pipeline_cmd = pipeline(cmd1, cmd2)

				output = read(run(pipeline_cmd), String)

				print(output)

				args_input = take!(input_channel)
				args = split(args_input)

				iX, iY, iZ, fX, fY, fZ = parse.(Float64, args)

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				write_sol_in_cuboid(iX, iY, iZ, fX, fY, fZ)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 5
				println("write_sol_in_multiple_cuboids(iX, iY, iZ, fX, fY, fZ, Nz) has been selected.\n")
				println("Arguments are required.\n")
				println("Please enter them as iX   iY   iZ   fX   fY   fZ   Nz.\n")
				println("As a reference, the dimensions of the NVT simulation box are")
				print("0.0   0.0   0.0")

				gro_files = glob("*.gro")
				cmd1 = `tail -q -n 1 $gro_files`
				cmd2 = `uniq`

				pipeline_cmd = pipeline(cmd1, cmd2)

				output = read(run(pipeline_cmd), String)

				print(output)

				args_input = take!(input_channel)
            args = split(args_input)

				iX, iY, iZ, fX, fY, fZ, Nz = parse.(Float64, args)

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				write_sol_in_multiple_cuboids(iX, iY, iZ, fX, fY, fZ, Nz)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 6
				println("write_dipole_moment_sol() has been selected, no arguments are needed.")
				println("Please be aware that the atom charges have already been established in")
				println("the code based on the MD parameters. If any adjustments are required,")
				println("please make them accordingly.")

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				write_dipole_moment_sol()
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 7
				println("process_sph_dipoles_sol() has been selected, no arguments are needed.")
				println("Please be aware that the atom charges have already been established in")
				println("the code based on the MD parameters. If any adjustments are required,")
				println("please make them accordingly.")

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				process_sph_dipoles_sol()
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 8
				println("write_sol_theta_in_multiple_cuboids(iX, iY, iZ, fX, fY, fZ, Nz) has been selected.\n")
				println("Arguments are required.\n")
				println("Please be aware that the atom charges have already been established in")
				println("the code based on the MD parameters. If any adjustments are required,")
				println("please make them accordingly.")
				println("Please enter them as iX   iY   iZ   fX   fY   fZ   Nz.\n")
				println("As a reference, the dimensions of the NVT simulation box are")
				print("0.0   0.0   0.0")

				gro_files = glob("*.gro")
				cmd1 = `tail -q -n 1 $gro_files`
				cmd2 = `uniq`

				pipeline_cmd = pipeline(cmd1, cmd2)

				output = read(run(pipeline_cmd), String)

				print(output)

				args_input = take!(input_channel)
				args = split(args_input)

				iX, iY, iZ, fX, fY, fZ, Nz = parse.(Float64, args)

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				write_sol_theta_in_multiple_cuboids(iX, iY, iZ, fX, fY, fZ, Nz)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 9
				println("rdf_sol(sol_molecule_name, dr, max_distance, iX, iY, iZ, fX, fY, fZ) has been selected.\n")
				println("Arguments are required.\n")
				println("Please enter them as sol_molecule_name   dr   max_distance   iX   iY   iZ   fX   fY   fZ.\n")
				println("As a reference, the dimensions of the NVT simulation box are")
				print("0.0   0.0   0.0")

				gro_files = glob("*.gro")
				cmd1 = `tail -q -n 1 $gro_files`
				cmd2 = `uniq`

				pipeline_cmd = pipeline(cmd1, cmd2)

				output = read(run(pipeline_cmd), String)

				print(output)

				args_input = take!(input_channel)
				args = split(args_input)

				sol_molecule_name = args[1]
				dr, max_distance, iX, iY, iZ, fX, fY, fZ = parse.(Float64, args[2:end])

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				rdf_sol(sol_molecule_name, dr, max_distance, iX, iY, iZ, fX, fY, fZ)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 10
				println("rdf_sol_all(dr, max_distance, iX, iY, iZ, fX, fY, fZ) has been selected.\n")
				println("Arguments are required.\n")
				println("Please enter them as dr   max_distance   iX   iY   iZ   fX   fY   fZ.\n")
				println("As a reference, the dimensions of the NVT simulation box are")
				print("0.0   0.0   0.0")

				gro_files = glob("*.gro")
				cmd1 = `tail -q -n 1 $gro_files`
				cmd2 = `uniq`

				pipeline_cmd = pipeline(cmd1, cmd2)

				output = read(run(pipeline_cmd), String)

				print(output)

				args_input = take!(input_channel)
				args = split(args_input)

				dr, max_distance, iX, iY, iZ, fX, fY, fZ = parse.(Float64, args[1:end])

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				rdf_sol_all(dr, max_distance, iX, iY, iZ, fX, fY, fZ)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == nothing || function_number > 10
				println("INVALID INPUT! Type help or exit.")
			end
		end
	end
end


