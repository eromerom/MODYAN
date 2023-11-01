
using Printf
using Glob

function read_gromacs_file(filename)
	time = 0.0
	molecule_names = String[]
	atom_names = String[]
	atom_coords = Tuple{Float64, Float64, Float64}[]
	box_dimensions = (0.0, 0.0, 0.0)

	open(filename, "r") do file
		line = "" # So that it is available outside the while loop
		while !eof(file)
			line = readline(file)
			if occursin("t=", line)
				time = parse(Float64, match(r"t=\s*([\d.]+)", line).captures[1])
			elseif (occursin("UN", line) || occursin("SOL", line)) && length(line) >= 44
				molecule_name = line[1:8]
				atom_name = line[9:20]
				x = parse(Float64, line[21:28])
				y = parse(Float64, line[29:36])
				z = parse(Float64, line[37:44])
				molecule_name = replace(molecule_name, " " => "")
				atom_name = replace(atom_name, " " => "")
				push!(molecule_names, molecule_name)
				push!(atom_names, atom_name)
				push!(atom_coords, (x, y, z))
			end
		end
		box_dimensions = Tuple(parse.(Float64, split(line)))
	end
	return time, molecule_names, atom_names, atom_coords, box_dimensions
end

function process_gromacs_files()

	file_list = glob("*.gro")

	results = Dict{Float64, 
				Tuple{Vector{String}, Vector{String}, 
				Vector{Tuple{Float64, Float64, Float64}}, 
				Tuple{Float64, Float64, Float64}}}()

	for filename in file_list
    	time, molecule_names, atom_names, atom_coords, box_dimensions = read_gromacs_file(filename)
    	results[time] = (molecule_names, atom_names, atom_coords, box_dimensions)
	end


	open("store_coordinates.log", "w") do file
	
		write(file, "Subroutine for processing GROMACS files.\nInitializing parsing and storing molecule and atom names with their corresponding coordinates.\n\n")

		for time in sort(collect(keys(results)))
			molecule_names, atom_names, atom_coords, box_dimensions = results[time]
			unique_molecules = unique(molecule_names)
			write(file, @sprintf("Time: %.5f\n", time))
			write(file, @sprintf("Total number of molecules: %d\n", length(unique_molecules)))
			write(file, @sprintf("Total number of atoms: %d\n", length(atom_names)))
			write(file, @sprintf("Box dimensions: %.5f %.5f %.5f\n\n", box_dimensions...))
		end

		write(file, "Success in storing GROMACS file data. This is the end of the log file.\n")

	end


	open("store_coordinates_data_saved.log", "w") do file
		for time in sort(collect(keys(results)))
			molecule_names, atom_names, atom_coords, box_dimensions = results[time]
			write(file, "Time: $time\n")
			write(file, "Format= Molecule name, atom name, coordinates (x, y, z)\n")
			for i in 1:length(molecule_names)
				write(file, "$(molecule_names[i]),$(atom_names[i]),$(atom_coords[i])\n")
			end
		end
	end
	return results
end

