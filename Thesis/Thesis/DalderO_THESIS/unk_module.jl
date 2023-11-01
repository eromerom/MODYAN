 
using LinearAlgebra
using Plots
using LaTeXStrings


function counting_unk()
	results = process_gromacs_files()

	open("count_unk.log", "w") do io
		println(io, "Reading UNK molecules...\n")
		for time in sort(collect(keys(results)))
			molecule_names, _, _, _ = results[time]
			unk_molecules = [molecule_name for molecule_name
			in molecule_names if molecule_name[end-2:end-1] == "UN"]
			unique_unk_molecules = Set(unk_molecules)
			unk_count = length(unique_unk_molecules)
			println(io, "Time: $time, number of UNK molecules: $unk_count")
		end
		println(io, "\nEnd of files, reading complete.")
	end
end



function write_unk()
   results = process_gromacs_files()

   open("write_unk.log", "w") do io
      for time in sort(collect(keys(results)))
         molecule_names, atom_names, atom_coords, _ = results[time]
			unk_indices = [i for i in 1:length(molecule_names) if molecule_names[i][end-2:end-1] == "UN"]
			write(io, "Time: $time\n")
			write(io, "Format= Molecule name, atom name, coordinates (x, y, z)\n")
			for i in unk_indices
				write(io, "$(molecule_names[i]),$(atom_names[i]),$(atom_coords[i])\n")
			end
      end
   end
end



function centre_of_mass_unk(results)
	centre_of_mass_unk_dict = Dict{Float64, Vector{Tuple{String, Tuple{Float64, Float64, Float64}}}}()

	for time in sort(collect(keys(results)))
		molecule_names, atom_names, atom_coords, box_dimensions = results[time]

		unk_indices = [i for i in 1:length(molecule_names) if molecule_names[i][end-2:end-1] == "UN"]
		unk_molecules = unique([molecule_names[i] for i in unk_indices])

		com_unk = []
		for unk_molecule in unk_molecules
			indices = [i for i in unk_indices if molecule_names[i] == unk_molecule]

			total_mass = 0.0
			com_x, com_y, com_z = (0.0, 0.0, 0.0)
			for i in indices
				atom_name = atom_names[i]
				(x, y, z) = atom_coords[i]

				mass = 0.0
				if atom_name[1] == 'O'
					mass = 15.9990
				elseif atom_name[1] == 'H'
					mass = 1.0080
				elseif atom_name[1] == 'C'
					mass = 12.0110
				end

				total_mass += mass
				com_x += mass * x
				com_y += mass * y
				com_z += mass * z
			end
			com_x /= total_mass
			com_y /= total_mass
			com_z /= total_mass

			push!(com_unk, (unk_molecule, (com_x, com_y, com_z)))
		end
		centre_of_mass_unk_dict[time] = com_unk
	end
	return centre_of_mass_unk_dict
end

function write_com_unk()
	results = process_gromacs_files()
	centre_of_mass_unk_dict = centre_of_mass_unk(results)

	open("com_unk.log", "w") do io
		for time in sort(collect(keys(centre_of_mass_unk_dict)))
			com_list = centre_of_mass_unk_dict[time]
			write(io, "Time: $time\n")
			write(io, "Number of UNK molecules in snapshot: $(length(com_list))\n" )
			write(io, "Format= UNK molecule, centre of mass coordinates (x, y, z)\n")
			for (molecule_name, com) in com_list
				write(io, "$molecule_name,($(com[1]),$(com[2]),$(com[3]))\n")
			end
		end
	end
end


function count_unk_in_cuboid(centre_of_mass_unk_dict, iX, iY, iZ, fX, fY, fZ)
	count_unk_in_cuboid_dict = Dict{Float64, Int64}()

	for time in sort(collect(keys(centre_of_mass_unk_dict)))
		com_list = centre_of_mass_unk_dict[time]
		count = 0
		for (_, com) in com_list
			x, y, z = com
			if iX <= x <= fX && iY <= y <= fY && iZ <= z <= fZ
				count += 1
			end
		end
		count_unk_in_cuboid_dict[time] = count
	end

	return count_unk_in_cuboid_dict
end

function write_unk_in_cuboid(iX, iY, iZ, fX, fY, fZ)

	results = process_gromacs_files()
	centre_of_mass_unk_dict = centre_of_mass_unk(results)
	count_unk_in_cuboid_dict = count_unk_in_cuboid(centre_of_mass_unk_dict, iX, iY, iZ, fX, fY, fZ)

	open("count_unk_in_cuboid.log", "w") do io
		println(io, "Counting UNK molecules in the cuboid...\n")
		println(io, "Initial point of the cuboid (iX, iY, iZ): ($iX, $iY, $iZ)")
		println(io, "Final point of the cuboid (fX, fY, fZ): ($fX, $fY, $fZ)\n")

		for time in sort(collect(keys(count_unk_in_cuboid_dict)))
			count = count_unk_in_cuboid_dict[time]
			println(io, "Time: $time, number of UNK molecules in cuboid: $count")
		end

		total_count = sum(values(count_unk_in_cuboid_dict))
		num_time_steps = length(count_unk_in_cuboid_dict)
		average_count = total_count / num_time_steps
		println(io, "\nAverage number of UNK molecules in cuboid over time: $average_count")

		println(io, "\nEnd of files, counting complete.")
	end
end


function write_unk_in_multiple_cuboids(iX, iY, iZ, fX, fY, fZ, Nz)

	results = process_gromacs_files()
	centre_of_mass_unk_dict = centre_of_mass_unk(results)
	bar_data_time = Dict()
	bar_data_avg = Dict()

	open("count_unk_in_multiple_cuboids.log", "w") do io
		println(io, "Counting UNK molecules in multiple cuboids...\n")
		println(io, "Initial point of the cuboids (iX, iY, iZ): ($iX, $iY, $iZ)")
		println(io, "Final point of the cuboids (fX, fY, fZ): ($fX, $fY, $fZ)\n")

		deltaZ = (fZ - iZ) / Nz

		for n in 1:Nz
			local_iZ = iZ + (n - 1) * deltaZ
			local_fZ = iZ + n * deltaZ
			count_unk_in_cuboid_dict = count_unk_in_cuboid(centre_of_mass_unk_dict,
			iX, iY, local_iZ,
			fX, fY, local_fZ)

			println(io, "Cuboid $n\n")
			println(io, "Initial point of the cuboid (iX, iY, iZ): ($iX, $iY, $local_iZ)")
			println(io, "Final point of the cuboid (fX, fY, fZ): ($fX, $fY, $local_fZ)\n")

			for time in sort(collect(keys(count_unk_in_cuboid_dict)))
				count = count_unk_in_cuboid_dict[time]
				println(io, "Time: $time, number of UNK molecules in cuboid $n: $count")

				p = bar([n], [count], title = "Time: $time", xlabel = "z (nm)", ylabel = "Number of Molecules", legend = false)
            savefig(p, "UNK_cuboid_$(n)_$(time).pdf")

            if haskey(bar_data_time, time)
               push!(bar_data_time[time], count)
            else
               bar_data_time[time] = [count]
            end
			end

			total_count = sum(values(count_unk_in_cuboid_dict))
			num_time_steps = length(count_unk_in_cuboid_dict)
			average_count = total_count / num_time_steps
			println(io, "\nAverage number of UNK molecules in cuboid $n over time: $average_count\n")

			bar_data_avg[n] = average_count
		end

		println(io, "End of files, counting complete.")
	end
	p_avg = bar(collect(keys(bar_data_avg)), collect(values(bar_data_avg)), title
	= "Average Number of UNK Molecules Over Time", xlabel = "z (nm)", ylabel = "Average Number of Molecules", legend = false)
	savefig(p_avg, "UNK_cuboid_average.pdf")
end


function dipole_moment_unk(results, qc00, qc01, qc02, qo03, qh04, qh05, qh06,
qh07, qc0a, qc0b, qc0c, qc0d, qh0e, qh0f, qh0g, qh0h, qh0i, qh0j)

	dipole_moment_unk_dict = Dict{Float64, Vector{Tuple{String, Tuple{Float64, Float64, Float64}}}}()

	for time in sort(collect(keys(results)))
		molecule_names, atom_names, atom_coords, box_dimensions = results[time]

		unk_indices = [i for i in 1:length(molecule_names) if molecule_names[i][end-2:end-1] == "UN"]
		unk_molecules = unique([molecule_names[i] for i in unk_indices])

		dipole_moments = []
		for unk_molecule in unk_molecules
			indices = [i for i in unk_indices if molecule_names[i] == unk_molecule]

			dip_x, dip_y, dip_z = (0.0, 0.0, 0.0)
			for i in indices
				atom_name = atom_names[i]
				(x, y, z) = atom_coords[i]

				charge = 0.0
				if atom_name[1:3] == "C00"
					charge = qc00
				elseif atom_name[1:3] == "C01"
					charge = qc01
				elseif atom_name[1:3] == "C02"
					charge = qc02
				elseif atom_name[1:3] == "O03"
					charge = qo03
				elseif atom_name[1:3] == "H04"
					charge = qh04
				elseif atom_name[1:3] == "H05"
					charge = qh05
				elseif atom_name[1:3] == "H06"
					charge = qh06
				elseif atom_name[1:3] == "H07"
					charge = qh07
				elseif atom_name[1:3] == "C0A"
					charge = qc0a
				elseif atom_name[1:3] == "C0B"
					charge = qc0b
				elseif atom_name[1:3] == "C0C"
					charge = qc0c
				elseif atom_name[1:3] == "C0D"
					charge = qc0d
				elseif atom_name[1:3] == "H0E"
					charge = qh0e
				elseif atom_name[1:3] == "H0F"
					charge = qh0f
				elseif atom_name[1:3] == "H0G"
					charge = qh0g
				elseif atom_name[1:3] == "H0H"
					charge = qh0h
				elseif atom_name[1:3] == "H0I"
					charge = qh0i
				elseif atom_name[1:3] == "H0J"
					charge = qh0j
				end 

				dip_x += charge * x
				dip_y += charge * y
				dip_z += charge * z
			end

			push!(dipole_moments, (unk_molecule, (dip_x, dip_y, dip_z)))
		end
		dipole_moment_unk_dict[time] = dipole_moments
	end
	return dipole_moment_unk_dict
end

function write_dipole_moment_unk()
	qc00 = -0.9618222787*0.2882 # Replace with your desired charge
	qc01 = -0.9618222787*0.1654 # Replace with your desired charge
	qc02 = 0.9618222787*0.3031 # Replace with your desired charge
	qo03 = -0.9618222787*0.3752 # Replace with your desired charge
	qh04 = 0.9618222787*0.0922 # Replace with your desired charge
	qh05 = 0.9618222787*0.1654 # Replace with your desired charge
	qh06 = 0.9618222787*0.1341 # Replace with your desired charge
	qh07 = 0.9618222787*0.1341 # Replace with your desired charge
	qc0a = -0.8286951127*0.1548 # Replace with your desired charge
	qc0b = -0.8286951127*0.1552 # Replace with your desired charge
	qc0c = -0.8286951127*0.2352 # Replace with your desired charge
	qc0d = -0.8286951127*0.2352 # Replace with your desired charge
	qh0e = 0.8286951127*0.1346 # Replace with your desired charge
	qh0f = 0.8286951127*0.1277 # Replace with your desired charge
	qh0g = 0.8286951127*0.1277 # Replace with your desired charge
	qh0h = 0.8286951127*0.1345 # Replace with your desired charge
	qh0i = 0.8286951127*0.1279 # Replace with your desired charge
	qh0j = 0.8286951127*0.1279 # Replace with your desired charge

	results = process_gromacs_files()
	dipole_moment_unk_dict = dipole_moment_unk(results, qc00, qc01, qc02, qo03, qh04, 
	qh05, qh06, qh07, qc0a, qc0b, qc0c, qc0d, qh0e, qh0f, qh0g, qh0h, qh0i, qh0j)

	open("dipole_moment_unk.log", "w") do io
 		for time in sort(collect(keys(dipole_moment_unk_dict)))
			dipole_list = dipole_moment_unk_dict[time]
			write(io, "Time: $time\n")
			write(io, "Number of UNK molecules in snapshot: $(length(dipole_list))\n")
			write(io, "Format= UNK molecule, dipole moment (x, y, z)\n")
			for (molecule_name, dipole) in dipole_list
				write(io, "$molecule_name,($(dipole[1]),$(dipole[2]),$(dipole[3]))\n")
			end
		end
	end
end


function cart_to_sph(cart_dipole::Tuple{Float64, Float64, Float64})
	x, y, z = cart_dipole
	r = norm(cart_dipole)
	theta = rad2deg(acos(z / r))
	phi = rad2deg(atan(y, x))
	phi = phi >= 0 ? phi : phi + 360
	return [r, theta, phi]
end

function store_sph_dipoles_unk(dipole_moment_unk_dict::Dict)
	sph_dipoles_unk = Dict{Float64, Dict{String, Vector{Float64}}}()

	for time in sort(collect(keys(dipole_moment_unk_dict)))
		dipole_list = dipole_moment_unk_dict[time]
		sph_dipoles_unk[time] = Dict{String, Vector{Float64}}()

		for (molecule_name, cart_dipole) in dipole_list
			sph_dipoles_unk[time][molecule_name] = cart_to_sph(cart_dipole)
		end
	end

	return sph_dipoles_unk
end

function write_sph_dipoles_unk(sph_dipoles_unk::Dict)
	open("sph_dipoles_unk.log", "w") do io
		for time in sort(collect(keys(sph_dipoles_unk)))
			sph_dipole_list = sph_dipoles_unk[time]
			write(io, "Time: $time\n")
			write(io, "Number of UNK molecules in snapshot: $(length(sph_dipole_list))\n")
			write(io, "Format= UNK molecule, spherical dipole moment (r, theta, phi)\n")

			sorted_molecules = sort(collect(keys(sph_dipole_list)))

			for molecule_name in sorted_molecules
				sph_dipole = sph_dipole_list[molecule_name]
				write(io, "$molecule_name,($(sph_dipole[1]),$(sph_dipole[2]),$(sph_dipole[3]))\n")
			end
		end
	end
end

function process_sph_dipoles_unk()
	qc00 = -0.9618222787*0.2882 # Replace with your desired charge
	qc01 = -0.9618222787*0.1654 # Replace with your desired charge
	qc02 = 0.9618222787*0.3031 # Replace with your desired charge
	qo03 = -0.9618222787*0.3752 # Replace with your desired charge
	qh04 = 0.9618222787*0.0922 # Replace with your desired charge
	qh05 = 0.9618222787*0.1654 # Replace with your desired charge
	qh06 = 0.9618222787*0.1341 # Replace with your desired charge
	qh07 = 0.9618222787*0.1341 # Replace with your desired charge
	qc0a = -0.8286951127*0.1548 # Replace with your desired charge
	qc0b = -0.8286951127*0.1552 # Replace with your desired charge
	qc0c = -0.8286951127*0.2352 # Replace with your desired charge
	qc0d = -0.8286951127*0.2352 # Replace with your desired charge
	qh0e = 0.8286951127*0.1346 # Replace with your desired charge
	qh0f = 0.8286951127*0.1277 # Replace with your desired charge
	qh0g = 0.8286951127*0.1277 # Replace with your desired charge
	qh0h = 0.8286951127*0.1345 # Replace with your desired charge
	qh0i = 0.8286951127*0.1279 # Replace with your desired charge
	qh0j = 0.8286951127*0.1279 # Replace with your desired charge

	results = process_gromacs_files()
	dipole_moment_unk_dict = dipole_moment_unk(results, qc00, qc01, qc02, qo03, qh04,
	qh05, qh06, qh07, qc0a, qc0b, qc0c, qc0d, qh0e, qh0f, qh0g, qh0h, qh0i, qh0j)
	sph_dipoles_unk = store_sph_dipoles_unk(dipole_moment_unk_dict)
	write_sph_dipoles_unk(sph_dipoles_unk)
end


function count_unk_in_cuboid_theta(centre_of_mass_unk_dict, sph_dipoles_unk,
iX, iY, iZ, fX, fY, fZ)

	count_unk_in_cuboid_theta_dict = Dict{Float64, Vector{Int}}()

	for time in sort(collect(keys(centre_of_mass_unk_dict)))
		com_list = centre_of_mass_unk_dict[time]
		theta_count = zeros(Int, 181)
		com_indices = []

		for (idx, com) in com_list
 			x, y, z = com
 			if iX <= x <= fX && iY <= y <= fY && iZ <= z <= fZ
 				push!(com_indices, idx)
			end
		end

 		sph_dipole_list = sph_dipoles_unk[time]

 		sorted_molecules = sort(collect(keys(sph_dipole_list)))

		for molecule_name in sorted_molecules
 			if molecule_name in com_indices
				sph_dipole = sph_dipole_list[molecule_name]
				theta = sph_dipole[2]
				theta_int = round(Int, theta)
				theta_count[theta_int + 1] += 1
			end
		end

		count_unk_in_cuboid_theta_dict[time] = theta_count
	end
	return count_unk_in_cuboid_theta_dict
end

function count_unk_in_cuboid_theta_dict(centre_of_mass_unk_dict, sph_dipoles_dict, iX, iY, iZ, fX, fY, fZ)
	qc00 = -0.9618222787*0.2882 # Replace with your desired charge
	qc01 = -0.9618222787*0.1654 # Replace with your desired charge
	qc02 = 0.9618222787*0.3031 # Replace with your desired charge
	qo03 = -0.9618222787*0.3752 # Replace with your desired charge
	qh04 = 0.9618222787*0.0922 # Replace with your desired charge
	qh05 = 0.9618222787*0.1654 # Replace with your desired charge
	qh06 = 0.9618222787*0.1341 # Replace with your desired charge
	qh07 = 0.9618222787*0.1341 # Replace with your desired charge
	qc0a = -0.8286951127*0.1548 # Replace with your desired charge
	qc0b = -0.8286951127*0.1552 # Replace with your desired charge
	qc0c = -0.8286951127*0.2352 # Replace with your desired charge
	qc0d = -0.8286951127*0.2352 # Replace with your desired charge
	qh0e = 0.8286951127*0.1346 # Replace with your desired charge
	qh0f = 0.8286951127*0.1277 # Replace with your desired charge
	qh0g = 0.8286951127*0.1277 # Replace with your desired charge
	qh0h = 0.8286951127*0.1345 # Replace with your desired charge
	qh0i = 0.8286951127*0.1279 # Replace with your desired charge
	qh0j = 0.8286951127*0.1279 # Replace with your desired charge

    theta_counts = Dict()
    for (time, unk_data) in centre_of_mass_unk_dict
	 	  dipole_moment_unk_dict = dipole_moment_unk(results,
		  qc00, qc01, qc02, qo03, qh04, qh05, qh06, qh07, qc0a,
		  qc0b, qc0c, qc0d, qh0e, qh0f, qh0g, qh0h, qh0i, qh0j)
        sph_dipoles_unk = store_sph_dipoles_unk(dipole_moment_unk_dict)
        theta_count = count_unk_in_cuboid_theta(unk_data, sph_dipoles_unk, iX, iY, iZ, fX, fY, fZ)
        theta_counts[time] = theta_count
    end
    return theta_counts
end

function write_unk_theta_in_multiple_cuboids(iX, iY, iZ, fX, fY, fZ, Nz)
	qc00 = -0.9618222787*0.2882 # Replace with your desired charge
	qc01 = -0.9618222787*0.1654 # Replace with your desired charge
	qc02 = 0.9618222787*0.3031 # Replace with your desired charge
	qo03 = -0.9618222787*0.3752 # Replace with your desired charge
	qh04 = 0.9618222787*0.0922 # Replace with your desired charge
	qh05 = 0.9618222787*0.1654 # Replace with your desired charge
	qh06 = 0.9618222787*0.1341 # Replace with your desired charge
	qh07 = 0.9618222787*0.1341 # Replace with your desired charge
	qc0a = -0.8286951127*0.1548 # Replace with your desired charge
	qc0b = -0.8286951127*0.1552 # Replace with your desired charge
	qc0c = -0.8286951127*0.2352 # Replace with your desired charge
	qc0d = -0.8286951127*0.2352 # Replace with your desired charge
	qh0e = 0.8286951127*0.1346 # Replace with your desired charge
	qh0f = 0.8286951127*0.1277 # Replace with your desired charge
	qh0g = 0.8286951127*0.1277 # Replace with your desired charge
	qh0h = 0.8286951127*0.1345 # Replace with your desired charge
	qh0i = 0.8286951127*0.1279 # Replace with your desired charge
	qh0j = 0.8286951127*0.1279 # Replace with your desired charge

    results = process_gromacs_files()
    centre_of_mass_unk_dict = centre_of_mass_unk(results)
	 dipole_moment_unk_dict = dipole_moment_unk(results, qc00, 
	 qc01, qc02, qo03, qh04, qh05, qh06, qh07, qc0a, qc0b, qc0c,
	 qc0d, qh0e, qh0f, qh0g, qh0h, qh0i, qh0j)
    sph_dipoles_dict = store_sph_dipoles_unk(dipole_moment_unk_dict)
	 bar_data_avg = Dict()

    open("count_unk_theta_in_multiple_cuboids.log", "w") do io
        println(io, "Counting UNK molecules in multiple cuboids with theta populations...\n")
        println(io, "Initial point of the cuboids (iX, iY, iZ): ($iX, $iY, $iZ)")
        println(io, "Final point of the cuboids (fX, fY, fZ): ($fX, $fY, $fZ)\n")

        deltaZ = (fZ - iZ) / Nz

        for n in 1:Nz
            local_iZ = iZ + (n - 1) * deltaZ
            local_fZ = iZ + n * deltaZ
            count_unk_in_cuboid_theta_dict = count_unk_in_cuboid_theta(centre_of_mass_unk_dict, sph_dipoles_dict,
            iX, iY, local_iZ,
            fX, fY, local_fZ)

				theta_count_avg = zeros(Float64, 181)

            println(io, "Cuboid $n\n")
            println(io, "Initial point of the cuboid (iX, iY, iZ): ($iX, $iY, $local_iZ)")
            println(io, "Final point of the cuboid (fX, fY, fZ): ($fX, $fY, $local_fZ)\n")

            for time in sort(collect(keys(count_unk_in_cuboid_theta_dict)))
                theta_counts = count_unk_in_cuboid_theta_dict[time]
					 theta_count_avg .+= theta_counts
                println(io, "Time: $time, theta counts in cuboid $n: $theta_counts")
            end

				theta_count_avg ./= length(count_unk_in_cuboid_theta_dict)
				bar_data_avg[n] = theta_count_avg

            total_count = sum(values(count_unk_in_cuboid_theta_dict))
            num_time_steps = length(count_unk_in_cuboid_theta_dict)
            average_count = total_count / num_time_steps
            println(io, "\nAverage number of UNK molecules in cuboid $n over time: $average_count\n")
        end

        println(io, "End of files, counting complete.")
    end
	 for (cuboid, theta_count_avg) in bar_data_avg
    	p_avg = bar(0:180, theta_count_avg, title = "Cuboid: $cuboid",
      xlabel = L"\theta", ylabel = "Average Count", legend = false)
      savefig(p_avg, "UNK_theta_cuboid_avg_$(cuboid).pdf")
    end
end


function distance_between(center_of_mass1::Tuple{Float64, Float64, Float64},
									center_of_mass2::Tuple{Float64, Float64, Float64})::Float64
	return norm(center_of_mass1 .- center_of_mass2)
end

function rdf_unk(unk_molecule_name, dr, max_distance, iX, iY, iZ, fX, fY, fZ)
	results = process_gromacs_files()
	centre_of_mass_dict = centre_of_mass_unk(results)

	Lx = fX - iX
	Ly = fY - iY
	Lz = fZ - iZ
	V_total = Lx * Ly * Lz

	time_step_plots = Dict{Int, Any}()
	time_step_rdfs = Dict{Int, Dict{Float64, Float64}}()

   total_rdf = Dict{Float64, Float64}()
   for r in 0:dr:max_distance
      total_rdf[r] = 0.0
   end

   translations = [(dx*Lx, dy*Ly, dz*Lz) for dx in -1:1 for dy in -1:1 for dz in -1:1]

	open("rdf_unk.log", "w") do io
		write(io, "Generating RDF plots...\n")
		for time in sort(collect(keys(centre_of_mass_dict)))

			rdf_histogram = Dict{Float64, Float64}()
			for r in 0:dr:max_distance
				rdf_histogram[r] = 0.0
			end

			com_list = centre_of_mass_dict[time]
			unk_centre_of_mass = [com for (molecule_name, com) 
			in com_list if molecule_name == unk_molecule_name][1]

			write(io, "Time: $time\n")

			for (molecule_name, com) in com_list
				if molecule_name == unk_molecule_name
					continue
				end

            for translation in translations
               translated_com = com .+ translation
               distance = distance_between(unk_centre_of_mass, translated_com)

					r_bin = floor(distance / dr) * dr
					if r_bin in keys(rdf_histogram)
						rdf_histogram[r_bin] += 1
					end
				end
			end

			time_step_rdf = copy(rdf_histogram)
			total_molecules = length(com_list)
			for r in keys(time_step_rdf)
				time_step_rdf[r] /= ((4/3) * pi * ((r + dr)^3 - r^3))
				time_step_rdf[r] /= (total_molecules / V_total)
				total_rdf[r] += time_step_rdf[r]
			end

			time_step_rdfs[time] = time_step_rdf

			x_values = collect(keys(time_step_rdf))
			y_values = collect(values(time_step_rdf))
			time_step_plots[time] = scatter(x_values, y_values, xlabel="Distance",
			ylabel="RDF", title="Radial Distribution Function (Time: $time)", legend=false)
			savefig(time_step_plots[time], "rdf_plot_$time.pdf")
		end
		write(io, "RDF plots generated, check for PDF files in current directory.")
	end

   num_snapshots = length(keys(centre_of_mass_dict))
   average_rdf = Dict{Float64, Float64}()
   for r in keys(total_rdf)
      average_rdf[r] = total_rdf[r] / num_snapshots
   end

   x_values = collect(keys(average_rdf))
   y_values = collect(values(average_rdf))
   average_rdf_plot = scatter(x_values, y_values, xlabel="Distance",
   ylabel="RDF", title="Average Radial Distribution Function", legend=false)
   savefig(average_rdf_plot, "average_rdf_plot.pdf")

	return time_step_rdfs
end


function rdf_unk_all(dr, max_distance, iX, iY, iZ, fX, fY, fZ)
	results = process_gromacs_files()
	centre_of_mass_dict = centre_of_mass_unk(results)

	Lx = fX - iX
	Ly = fY - iY
	Lz = fZ - iZ
	V_total = Lx * Ly * Lz

	time_step_plots = Dict{Int, Any}()
	time_step_rdfs = Dict{Int, Dict{Float64, Float64}}()

	total_rdf = Dict{Float64, Float64}()
	for r in 0:dr:max_distance
		total_rdf[r] = 0.0
	end

	translations = [(dx*Lx, dy*Ly, dz*Lz) for dx in -1:1 for dy in -1:1 for dz in -1:1]

	open("rdf_unk_all.log", "w") do io
		write(io, "Generating RDF plots...\n")
		for time in sort(collect(keys(centre_of_mass_dict)))

			rdf_histogram = Dict{Float64, Float64}()
			for r in 0:dr:max_distance
				rdf_histogram[r] = 0.0
			end

			com_list = centre_of_mass_dict[time]

			write(io, "Time: $time\n")

			for (unk_molecule_name, unk_centre_of_mass) in com_list

				for (molecule_name, com) in com_list
					if molecule_name == unk_molecule_name
						continue
					end

					for translation in translations
						translated_com = com .+ translation
						distance = distance_between(unk_centre_of_mass, translated_com)

						r_bin = floor(distance / dr) * dr
						if r_bin in keys(rdf_histogram)
							rdf_histogram[r_bin] += 1
						end
					end
				end
			end

			time_step_rdf = copy(rdf_histogram)
			total_molecules = length(com_list)
			for r in keys(time_step_rdf)
				time_step_rdf[r] /= (2 * total_molecules)
				time_step_rdf[r] /= ((4/3) * pi * ((r + dr)^3 - r^3))
				time_step_rdf[r] /= (total_molecules / V_total)
				time_step_rdfs[time] = time_step_rdf
			end

			x_values = collect(keys(time_step_rdfs[time]))
			y_values = collect(values(time_step_rdfs[time]))
			time_step_plots[time] = scatter(x_values, y_values, xlabel="Distance",
			ylabel="RDF", title="Radial Distribution Function (Time: $time)", legend=false)
			savefig(time_step_plots[time], "rdf_plot_all_$time.pdf")
		end
		write(io, "RDF plots generated, check for PDF files in current directory.")
	end

	num_snapshots = length(keys(centre_of_mass_dict))
	average_rdf = Dict{Float64, Float64}()
	for r in keys(total_rdf)
		for time in keys(time_step_rdfs)
			total_rdf[r] += time_step_rdfs[time][r]
		end
		average_rdf[r] = total_rdf[r] / num_snapshots
	end

	x_values = collect(keys(average_rdf))
	y_values = collect(values(average_rdf))
	average_rdf_plot = scatter(x_values, y_values, xlabel="Distance",
	ylabel="RDF", title="Average Radial Distribution Function", legend=false)
	savefig(average_rdf_plot, "average_rdf_plot_all.pdf")

	return time_step_rdfs
end

