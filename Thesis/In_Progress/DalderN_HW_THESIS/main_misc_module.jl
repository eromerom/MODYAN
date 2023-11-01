
function show_help_misc()
	println("""
	┌─────────────────────────────────────────────────────────────────────────────────────┐
	│                             List of Available Functions                             │
	├─────────────────────────────────────────────────────────────────────────────────────┤
	│ 1: Generate RDF plots for UNK-SOL (single mulecule).                                │
	│ 2: Generate a list of UNK molecules solvated within a distance by N water molecules.│
	│ 3: Electric Field in a point between two atoms over time.                           │
	│ 4: Generate RDF plots for UNK-SOL (all molecules).                                  │
	│ 5: Slice analyzer of electric fields.                                               │
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

function main_misc()

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
			show_help_misc()
		else

			function_number = tryparse(Int, input)

			if function_number == 1
				println("rdf_misc(unk_molecule_name, dr, max_distance, iX, iY, iZ, fX, fY, fZ) has been selected.\n")
				println("Arguments are required.\n")
				println("Please enter them as unk_molecule_name   dr   max_distance   iX   iY   iZ   fX   fY   fZ.\n")
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

				unk_molecule_name = args[1]
				dr, max_distance, iX, iY, iZ, fX, fY, fZ = parse.(Float64, args[2:end])

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				rdf_misc(unk_molecule_name, dr, max_distance, iX, iY, iZ, fX, fY, fZ)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 2
				println("solvation_searcher(max_distance, NW) has been selected.\n")
				println("Arguments are required.\n")
				println("Please enter them as max_distance   NW.\n")
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

				max_distance, NW = parse.(Float64, args)

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				solvation_searcher(max_distance, NW)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 3
				println("electric_field_at_point(molecule_name, atom_1, atom_2, iX, iY, iZ, fX, fY, fZ) has been selected.\n")
				println("Arguments are required.\n")
				println("Please enter them as molecule_name   atom_1   atom_2   iX   iY   iZ   fX   fY   fZ.\n")
				println("As a reference, atoms corresponding to the C-O bond are")

            cmd = `echo "O04   C07"`

            output = read(run(cmd), String)

            print(output)

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

				molecule_name = args[1]
				atom_1 = args[2]
				atom_2 = args[3]
				iX, iY, iZ, fX, fY, fZ = parse.(Float64, args[4:end])

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				electric_field_at_point(molecule_name, atom_1, atom_2, iX, iY, iZ, fX, fY, fZ)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 4
				println("rdf_misc_all(dr, max_distance, iX, iY, iZ, fX, fY, fZ) has been selected.\n")
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
				rdf_misc_all(dr, max_distance, iX, iY, iZ, fX, fY, fZ)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == 5
				println("write_angle_in_multiple_slices(molecule_type, atom_1, atom_2, atom_3, iX, iY, iZ, fX, fY, fZ, Nz) has been selected.\n")
				println("Arguments are required.\n")
				println("Please enter them as molecule_type   atom_1   atom_2  atom_3  iX   iY   iZ   fX   fY   fZ   Nz.\n")
				println("As a reference, the dimensions of the NVT simulation box are")
				print("0.0   0.0   0.0")

				gro_files = glob("*.gro")
				cmd1 = `tail -q -n 1 $gro_files`
				cmd2 = `uniq`

				pipeline_cmd = pipeline(cmd1, cmd2)

				output = read(run(pipeline_cmd), String)

				print(output)

            println("As a reference, atoms corresponding to the C-O bond are")

            cmd = `echo "O04   C07"`

            output = read(run(cmd), String)

            print(output)

				args_input = take!(input_channel)
				args = split(args_input)

				molecule_type = args[1]
				atom_1 = args[2]
				atom_2 = args[3]
				atom_3 = args[4]
				iX, iY, iZ, fX, fY, fZ, Nz = parse.(Float64, args[5:end])

				stop_channel = Channel(1)
				progress_task = @async show_progress(stop_channel)
				write_angle_in_multiple_slices(molecule_type, atom_1, atom_2, atom_3, iX, iY, iZ, fX, fY, fZ, Nz)
				put!(stop_channel, true)
				@async progress_task

			elseif function_number == nothing || function_number > 5
				println("INVALID INPUT! Type help or exit.")
			end
		end
	end
end


