
using LinearAlgebra
using Plots
using LaTeXStrings


function counting_sol()
	results = process_gromacs_files()

	open("count_sol.log", "w") do io
		println(io, "Reading SOL molecules...\n")
		for time in sort(collect(keys(results)))
			molecule_names, _, _, _ = results[time]
			sol_molecules = [molecule_name for molecule_name
			in molecule_names if molecule_name[end-2:end] == "SOL"]
			unique_sol_molecules = Set(sol_molecules)
			sol_count = length(unique_sol_molecules)
			println(io, "Time: $time, number of SOL molecules: $sol_count")
		end
		println(io, "\nEnd of files, reading complete.")
	end
end



function write_sol()
   results = process_gromacs_files()

   open("write_sol.log", "w") do io
      for time in sort(collect(keys(results)))
         molecule_names, atom_names, atom_coords, box_dimensions = results[time]
			sol_indices = [i for i in 1:length(molecule_names) if molecule_names[i][end-2:end] == "SOL"]
			write(io, "Time: $time\n")
			write(io, "Format= Molecule name, atom name, coordinates (x, y, z)\n")
			for i in sol_indices
				write(io, "$(molecule_names[i]),$(atom_names[i]),$(atom_coords[i])\n")
			end
      end
   end
end



function centre_of_mass_sol(results)
	centre_of_mass_sol_dict = Dict{Float64, Vector{Tuple{String, Tuple{Float64, Float64, Float64}}}}()

	for time in sort(collect(keys(results)))
		molecule_names, atom_names, atom_coords, box_dimensions = results[time]

		sol_indices = [i for i in 1:length(molecule_names) if molecule_names[i][end-2:end] == "SOL"]
		sol_molecules = unique([molecule_names[i] for i in sol_indices])

		com_sol = []
		for sol_molecule in sol_molecules
			indices = [i for i in sol_indices if molecule_names[i] == sol_molecule]

			total_mass = 0.0
			com_x, com_y, com_z = (0.0, 0.0, 0.0)
			for i in indices
				atom_name = atom_names[i]
				(x, y, z) = atom_coords[i]

				mass = 0.0
				if atom_name[1:2] == "OW"
					mass = 15.99940
				elseif atom_name[1:2] == "HW"
					mass = 1.008
				end

				total_mass += mass
				com_x += mass * x
				com_y += mass * y
				com_z += mass * z
			end
			com_x /= total_mass
			com_y /= total_mass
			com_z /= total_mass

			push!(com_sol, (sol_molecule, (com_x, com_y, com_z)))
		end
		centre_of_mass_sol_dict[time] = com_sol
	end
	return centre_of_mass_sol_dict
end

function write_com_sol()
	results = process_gromacs_files()
	centre_of_mass_sol_dict = centre_of_mass_sol(results)

	open("com_sol.log", "w") do io
		for time in sort(collect(keys(centre_of_mass_sol_dict)))
			com_list = centre_of_mass_sol_dict[time]
			write(io, "Time: $time\n")
			write(io, "Number of SOL molecules in snapshot: $(length(com_list))\n" )
			write(io, "Format= SOL molecule, centre of mass coordinates (x, y, z)\n")
			for (molecule_name, com) in com_list
				write(io, "$molecule_name,($(com[1]),$(com[2]),$(com[3]))\n")
			end
		end
	end
end


function count_sol_in_cuboid(centre_of_mass_sol_dict, iX, iY, iZ, fX, fY, fZ)
	count_sol_in_cuboid_dict = Dict{Float64, Int64}()

	for time in sort(collect(keys(centre_of_mass_sol_dict)))
		com_list = centre_of_mass_sol_dict[time]
		count = 0
		for (_, com) in com_list
			x, y, z = com
			if iX <= x <= fX && iY <= y <= fY && iZ <= z <= fZ
				count += 1
			end
		end
		count_sol_in_cuboid_dict[time] = count
	end

	return count_sol_in_cuboid_dict
end

function write_sol_in_cuboid(iX, iY, iZ, fX, fY, fZ)

	results = process_gromacs_files()
	centre_of_mass_sol_dict = centre_of_mass_sol(results)
	count_sol_in_cuboid_dict = count_sol_in_cuboid(centre_of_mass_sol_dict, iX, iY, iZ, fX, fY, fZ)

	open("count_sol_in_cuboid.log", "w") do io
		println(io, "Counting SOL molecules in the cuboid...\n")
		println(io, "Initial point of the cuboid (iX, iY, iZ): ($iX, $iY, $iZ)")
		println(io, "Final point of the cuboid (fX, fY, fZ): ($fX, $fY, $fZ)\n")

		for time in sort(collect(keys(count_sol_in_cuboid_dict)))
			count = count_sol_in_cuboid_dict[time]
			println(io, "Time: $time, number of SOL molecules in cuboid: $count")
		end

		total_count = sum(values(count_sol_in_cuboid_dict))
		num_time_steps = length(count_sol_in_cuboid_dict)
		average_count = total_count / num_time_steps
		println(io, "\nAverage number of SOL molecules in cuboid over time: $average_count")

		println(io, "\nEnd of files, counting complete.")
	end
end


function write_sol_in_multiple_cuboids(iX, iY, iZ, fX, fY, fZ, Nz)

	results = process_gromacs_files()
	centre_of_mass_sol_dict = centre_of_mass_sol(results)
	bar_data_time = Dict()
	bar_data_avg = Dict()

	open("count_sol_in_multiple_cuboids.log", "w") do io
		println(io, "Counting SOL molecules in multiple cuboids...\n")
		println(io, "Initial point of the cuboids (iX, iY, iZ): ($iX, $iY, $iZ)")
		println(io, "Final point of the cuboids (fX, fY, fZ): ($fX, $fY, $fZ)\n")

		deltaZ = (fZ - iZ) / Nz

		for n in 1:Nz
			local_iZ = iZ + (n - 1) * deltaZ
			local_fZ = iZ + n * deltaZ
			count_sol_in_cuboid_dict = count_sol_in_cuboid(centre_of_mass_sol_dict,
			iX, iY, local_iZ,
			fX, fY, local_fZ)

			println(io, "Cuboid $n\n")
			println(io, "Initial point of the cuboid (iX, iY, iZ): ($iX, $iY, $local_iZ)")
			println(io, "Final point of the cuboid (fX, fY, fZ): ($fX, $fY, $local_fZ)\n")

			for time in sort(collect(keys(count_sol_in_cuboid_dict)))
				count = count_sol_in_cuboid_dict[time]
				println(io, "Time: $time, number of SOL molecules in cuboid $n: $count")

				p = bar([n], [count], title = "Time: $time", xlabel = "z (nm)", ylabel = "Number of Molecules", legend = false)
				savefig(p, "SOL_cuboid_$(n)_$(time).pdf")

				if haskey(bar_data_time, time)
					push!(bar_data_time[time], count)
				else
					bar_data_time[time] = [count]
				end
			end

			total_count = sum(values(count_sol_in_cuboid_dict))
			num_time_steps = length(count_sol_in_cuboid_dict)
			average_count = total_count / num_time_steps
			println(io, "\nAverage number of SOL molecules in cuboid $n over time: $average_count\n")

			bar_data_avg[n] = average_count
		end

		println(io, "End of files, counting complete.")
	end
	p_avg = bar(collect(keys(bar_data_avg)), collect(values(bar_data_avg)), title 
	= "", xlabel = "z (nm)", ylabel = "Average Number of Molecules", legend = false, ylim = (0, 11), xlim = (0, 9.5), color = :red )
	savefig(p_avg, "SOL_cuboid_average.pdf")
end


function dipole_moment_sol(results, qow, qhw)
	dipole_moment_sol_dict = Dict{Float64, Vector{Tuple{String, Tuple{Float64, Float64, Float64}}}}()

	for time in sort(collect(keys(results)))
		molecule_names, atom_names, atom_coords, box_dimensions = results[time]

		sol_indices = [i for i in 1:length(molecule_names) if molecule_names[i][end-2:end] == "SOL"]
		sol_molecules = unique([molecule_names[i] for i in sol_indices])

		dipole_moments = []
		for sol_molecule in sol_molecules
			indices = [i for i in sol_indices if molecule_names[i] == sol_molecule]

			dip_x, dip_y, dip_z = (0.0, 0.0, 0.0)
			for i in indices
				atom_name = atom_names[i]
				(x, y, z) = atom_coords[i]

				charge = 0.0
				if atom_name[1:2] == "OW"
					charge = qow
				elseif atom_name[1:2] == "HW"
					charge = qhw
				end

				dip_x += charge * x
				dip_y += charge * y
				dip_z += charge * z
			end

			push!(dipole_moments, (sol_molecule, (dip_x, dip_y, dip_z)))
		end
		dipole_moment_sol_dict[time] = dipole_moments
	end
	return dipole_moment_sol_dict
end

function write_dipole_moment_sol()
	qow = -0.7825325797*0.834 # Replace with your desired charge
	qhw = 0.7825325797*0.417  # Replace with your desired charge

	results = process_gromacs_files()
	dipole_moment_sol_dict = dipole_moment_sol(results, qow, qhw)

	open("dipole_moment_sol.log", "w") do io
 		for time in sort(collect(keys(dipole_moment_sol_dict)))
			dipole_list = dipole_moment_sol_dict[time]
			write(io, "Time: $time\n")
			write(io, "Number of SOL molecules in snapshot: $(length(dipole_list))\n")
			write(io, "Format= SOL molecule, dipole moment (x, y, z)\n")
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

function store_sph_dipoles_sol(dipole_moment_sol_dict::Dict)
	sph_dipoles_sol = Dict{Float64, Dict{String, Vector{Float64}}}()

	for time in sort(collect(keys(dipole_moment_sol_dict)))
		dipole_list = dipole_moment_sol_dict[time]
		sph_dipoles_sol[time] = Dict{String, Vector{Float64}}()

		for (molecule_name, cart_dipole) in dipole_list
			sph_dipoles_sol[time][molecule_name] = cart_to_sph(cart_dipole)
		end
	end

	return sph_dipoles_sol
end

function write_sph_dipoles_sol(sph_dipoles_sol::Dict)
	open("sph_dipoles_sol.log", "w") do io
		for time in sort(collect(keys(sph_dipoles_sol)))
			sph_dipole_list = sph_dipoles_sol[time]
			write(io, "Time: $time\n")
			write(io, "Number of SOL molecules in snapshot: $(length(sph_dipole_list))\n")
			write(io, "Format= SOL molecule, spherical dipole moment (r, theta, phi)\n")

			sorted_molecules = sort(collect(keys(sph_dipole_list)))

			for molecule_name in sorted_molecules
				sph_dipole = sph_dipole_list[molecule_name]
				write(io, "$molecule_name,($(sph_dipole[1]),$(sph_dipole[2]),$(sph_dipole[3]))\n")
			end
		end
	end
end

function process_sph_dipoles_sol()
	qow = -0.7825325797*0.834 # Replace with your desired charge
	qhw = 0.7825325797*0.417  # Replace with your desired charge

	results = process_gromacs_files()
	dipole_moment_sol_dict = dipole_moment_sol(results, qow, qhw)
	sph_dipoles_sol = store_sph_dipoles_sol(dipole_moment_sol_dict)
	write_sph_dipoles_sol(sph_dipoles_sol)
end


function count_sol_in_cuboid_theta(centre_of_mass_sol_dict, sph_dipoles_sol,
iX, iY, iZ, fX, fY, fZ)

	count_sol_in_cuboid_theta_dict = Dict{Float64, Vector{Int}}()

	for time in sort(collect(keys(centre_of_mass_sol_dict)))
		com_list = centre_of_mass_sol_dict[time]
		theta_count = zeros(Int, 181)
		com_indices = []

		for (idx, com) in com_list
 			x, y, z = com
 			if iX <= x <= fX && iY <= y <= fY && iZ <= z <= fZ
 				push!(com_indices, idx)
			end
		end

 		sph_dipole_list = sph_dipoles_sol[time]

 		sorted_molecules = sort(collect(keys(sph_dipole_list)))

		for molecule_name in sorted_molecules
 			if molecule_name in com_indices
				sph_dipole_sol = sph_dipole_list[molecule_name]
				theta = sph_dipole_sol[2]
				theta_int = round(Int, theta)
				theta_count[theta_int + 1] += 1
			end
		end
      	
		count_sol_in_cuboid_theta_dict[time] = theta_count
	end
	return count_sol_in_cuboid_theta_dict
end

function count_sol_in_cuboid_theta_dict(centre_of_mass_sol_dict, sph_dipoles_dict, iX, iY, iZ, fX, fY, fZ)
	qow = -0.7825325797*0.834 # Replace with your desired charge
	qhw = 0.7825325797*0.417  # Replace with your desired charge

    theta_counts = Dict()
    for (time, sol_data) in centre_of_mass_sol_dict
	 	  dipole_moment_sol_dict = dipole_moment_sol(results, qow, qhw)	
        sph_dipoles_sol = store_sph_dipoles_sol(dipole_moment_sol_dict)
        theta_count = count_sol_in_cuboid_theta(sol_data, sph_dipoles_sol, iX, iY, iZ, fX, fY, fZ)
        theta_counts[time] = theta_count
    end
    return theta_counts
end

function write_sol_theta_in_multiple_cuboids(iX, iY, iZ, fX, fY, fZ, Nz)
	qow = -0.7825325797*0.834 # Replace with your desired charge
	qhw = 0.7825325797*0.417  # Replace with your desired charge

    results = process_gromacs_files()
    centre_of_mass_sol_dict = centre_of_mass_sol(results)
	 dipole_moment_sol_dict = dipole_moment_sol(results, qow, qhw)
    sph_dipoles_dict = store_sph_dipoles_sol(dipole_moment_sol_dict)
	 bar_data_avg = Dict()

    open("count_sol_theta_in_multiple_cuboids.log", "w") do io
        println(io, "Counting SOL molecules in multiple cuboids with theta populations...\n")
        println(io, "Initial point of the cuboids (iX, iY, iZ): ($iX, $iY, $iZ)")
        println(io, "Final point of the cuboids (fX, fY, fZ): ($fX, $fY, $fZ)\n")

        deltaZ = (fZ - iZ) / Nz

        for n in 1:Nz
            local_iZ = iZ + (n - 1) * deltaZ
            local_fZ = iZ + n * deltaZ
            count_sol_in_cuboid_theta_dict = count_sol_in_cuboid_theta(centre_of_mass_sol_dict, sph_dipoles_dict,
            iX, iY, local_iZ,
            fX, fY, local_fZ)

				theta_count_avg = zeros(Float64, 181)

            println(io, "Cuboid $n\n")
            println(io, "Initial point of the cuboid (iX, iY, iZ): ($iX, $iY, $local_iZ)")
            println(io, "Final point of the cuboid (fX, fY, fZ): ($fX, $fY, $local_fZ)\n")

            for time in sort(collect(keys(count_sol_in_cuboid_theta_dict)))
                theta_counts = count_sol_in_cuboid_theta_dict[time]
					 theta_count_avg .+= theta_counts
                println(io, "Time: $time, theta counts in cuboid $n: $theta_counts")
            end

				theta_count_avg ./= length(count_sol_in_cuboid_theta_dict)
				bar_data_avg[n] = theta_count_avg

            total_count = sum(values(count_sol_in_cuboid_theta_dict))
            num_time_steps = length(count_sol_in_cuboid_theta_dict)
            average_count = total_count / num_time_steps
            println(io, "\nAverage number of SOL molecules in cuboid $n over time: $average_count\n")
        end

        println(io, "End of files, counting complete.")
    end
	 for (cuboid, theta_count_avg) in bar_data_avg
	 	p_avg = bar(0:180, theta_count_avg, title = "Cuboid: $cuboid", 
		xlabel = L"\theta", ylabel = "Average Count", legend = false)
		savefig(p_avg, "SOL_theta_cuboid_avg_$(cuboid).pdf")
	 end
end


function distance_between(center_of_mass1::Tuple{Float64, Float64, Float64},
									center_of_mass2::Tuple{Float64, Float64, Float64})::Float64
	return norm(center_of_mass1 .- center_of_mass2)
end

function rdf_sol(sol_molecule_name, dr, max_distance, iX, iY, iZ, fX, fY, fZ)
	results = process_gromacs_files()
	centre_of_mass_dict = centre_of_mass_sol(results)

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

	open("rdf_sol.log", "w") do io
		write(io, "Generating RDF plots...\n")
		for time in sort(collect(keys(centre_of_mass_dict)))

			rdf_histogram = Dict{Float64, Float64}()
			for r in 0:dr:max_distance
				rdf_histogram[r] = 0.0
			end

			com_list = centre_of_mass_dict[time]
			sol_centre_of_mass = [com for (molecule_name, com) 
			in com_list if molecule_name == sol_molecule_name][1]

			write(io, "Time: $time\n")

			for (molecule_name, com) in com_list
				if molecule_name == sol_molecule_name
					continue
				end

				for translation in translations
					translated_com = com .+ translation
					distance = distance_between(sol_centre_of_mass, translated_com)

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


function rdf_sol_all(dr, max_distance, iX, iY, iZ, fX, fY, fZ)
	results = process_gromacs_files()
	centre_of_mass_dict = centre_of_mass_sol(results)

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

	open("rdf_sol_all.log", "w") do io
		write(io, "Generating RDF plots...\n")
		for time in sort(collect(keys(centre_of_mass_dict)))

			rdf_histogram = Dict{Float64, Float64}()
			for r in 0:dr:max_distance
				rdf_histogram[r] = 0.0
			end

			com_list = centre_of_mass_dict[time]

			write(io, "Time: $time\n")

			for (sol_molecule_name, sol_centre_of_mass) in com_list

				for (molecule_name, com) in com_list
					if molecule_name == sol_molecule_name
						continue
					end

					for translation in translations
						translated_com = com .+ translation
						distance = distance_between(sol_centre_of_mass, translated_com)

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
