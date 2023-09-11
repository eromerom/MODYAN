 
using LinearAlgebra
using Plots
using LaTeXStrings
using Statistics
using KernelDensity
using CSV
using DataFrames

function cross_product_tuple(a::Tuple{Float64, Float64, Float64}, b::Tuple{Float64, Float64, Float64})
    x = a[2]*b[3] - a[3]*b[2]
    y = a[3]*b[1] - a[1]*b[3]
    z = a[1]*b[2] - a[2]*b[1]
    return (x, y, z)
end

function distance_between(center_of_mass1::Tuple{Float64, Float64, Float64},
									center_of_mass2::Tuple{Float64, Float64, Float64})::Float64
	return norm(center_of_mass1 .- center_of_mass2)
end

function rdf_misc(unk_molecule_name, dr, max_distance, iX, iY, iZ, fX, fY, fZ)
	results = process_gromacs_files()
	centre_of_mass_dict_unk = centre_of_mass_unk(results)

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

	open("rdf_misc.log", "w") do io
		write(io, "Generating RDF plots...\n")
		for time in sort(collect(keys(centre_of_mass_dict_unk)))

			rdf_histogram = Dict{Float64, Float64}()
			for r in 0:dr:max_distance
				rdf_histogram[r] = 0.0
			end

			com_list_unk = centre_of_mass_dict_unk[time]
			unk_centre_of_mass = [com for (molecule_name, com) 
			in com_list_unk if molecule_name == unk_molecule_name][1]

			write(io, "Time: $time\n")

			centre_of_mass_dict_sol = centre_of_mass_sol(results)
			com_list_sol = centre_of_mass_dict_sol[time]

			for (molecule_name, com) in com_list_sol
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
			total_molecules = length(com_list_sol)
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

   num_snapshots = length(keys(centre_of_mass_dict_unk))
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


function solvation_searcher(max_distance, NW)
	results = process_gromacs_files()
	centre_of_mass_dict_unk = centre_of_mass_unk(results)
	centre_of_mass_dict_sol = centre_of_mass_sol(results)
	found_molecules = false

	open("solvated_UNK.log", "w") do io
		for time in sort(collect(keys(centre_of_mass_dict_unk)))
			com_list_unk = centre_of_mass_dict_unk[time]
			com_list_sol = centre_of_mass_dict_sol[time]
			write(io, "Time: $time\n")

			for (molecule_name, com_unk) in com_list_unk
				nearby_sol_count = 0

				for (sol_molecule_name, com_sol) in com_list_sol
					distance = distance_between(com_unk, com_sol)

					if distance <= max_distance
						nearby_sol_count += 1
					end
				end

				if nearby_sol_count == NW
					write(io, "Molecule name: $molecule_name\n")
					found_molecules = true
				end
			end
		end

		if !found_molecules
			write(io, "No molecules with $NW neighbours within r=$max_distance found.")
		end
	end
end


function get_atomic_charges()
	return Dict(
		"OW" => -0.7825325797*0.834,
		"HW" => 0.7825325797*0.417,
		"C00" => -0.8879014056*0.1983,
		"C01" => -0.8879014056*0.1724,
		"C02" => 0.8879014056*0.1328,
		"H03" => 0.8879014056*0.1064,
		"N04" => -0.8879014056*0.6161,
		"H05" => 0.8879014056*0.1414,
		"H06" => 0.8879014056*0.1393,
		"H07" => 0.8879014056*0.1393,
		"H08" => 0.8879014056*0.3276,
		"C0A" => -1.08683116*0.1576,
		"C0B" => -1.08683116*0.1561,
		"C0C" => -1.08683116*0.2363,
		"C0D" => -1.08683116*0.2356,
		"H0E" => 1.08683116*0.1356,
		"H0F" => 1.08683116*0.1282,
		"H0G" => 1.08683116*0.1282,
		"H0H" => 1.08683116*0.1370,
		"H0I" => 1.08683116*0.1283,
		"H0J" => 1.08683116*0.1283,
							)
end

function middle_of_atoms(results, molecule_name, atom_1, atom_2, time)
	molecule_names, atom_names, atom_coords, box_dimensions = results[time]

	molecule_indices = findall(==(molecule_name), molecule_names)
	atom_1_index = findfirst(x -> startswith(x, atom_1), atom_names[molecule_indices])
	atom_2_index = findfirst(x -> startswith(x, atom_2), atom_names[molecule_indices])

	if atom_1_index !== nothing && atom_2_index !== nothing
		atom_1_coords = atom_coords[molecule_indices[atom_1_index]]
 		atom_2_coords = atom_coords[molecule_indices[atom_2_index]]

		middle_point = 0.5 .* (atom_1_coords .+ atom_2_coords)

		return middle_point
	end
end

function electric_field_at_point(molecule_name, atom_1, atom_2, atom_3, iX, iY, iZ, fX, fY, fZ)
	k = 1.43996454 # Coulomb's constant in (eV nm^2 {elemental charge}^(-2))
   Lx = fX - iX
   Ly = fY - iY
   Lz = fZ - iZ

	atomic_charges = get_atomic_charges()

	electric_field_results = Dict{Float64, Tuple{Float64, Float64, Float64}}()
	angle_results = Dict{Float64, Float64}()
	projection_results = Dict{Float64, Tuple{Float64, Float64, Float64}}()
	magnitude_projection_results = Dict{Float64, Float64}()
	translations = [(dx*Lx, dy*Ly, dz*Lz) for dx in -1:1 for dy in -1:1 for dz in -1:1]
	results = process_gromacs_files()

	for time in sort(collect(keys(results)))
		middle_point = middle_of_atoms(results, molecule_name, atom_1, atom_2, time)

		Electric_Field = (0.0, 0.0, 0.0)
		molecule_names, atom_names, atom_coords, box_dimensions = results[time]

		target_atom_names = atom_names
		target_atom_coords = atom_coords
		target_molecule_names = molecule_names

		for (molecule, atom_name, (x, y, z)) in zip(target_molecule_names, target_atom_names, target_atom_coords)

			if molecule == molecule_name continue end

			q = nothing

			for key in keys(atomic_charges)
				if startswith(atom_name, key)
					q = atomic_charges[key]
					break
				end
			end

			if q === nothing
				continue
			end
			
			for translation in translations
				translated_atom_coords = (x, y, z) .+ translation
				r = translated_atom_coords .- middle_point
				d = norm(r)
				E_i = (k * q) .* r ./ d^3
				Electric_Field = Electric_Field .+ E_i
			end
		end
		
		atom_1_index = findfirst(x -> startswith(x, atom_1), target_atom_names)
		atom_2_index = findfirst(x -> startswith(x, atom_2), target_atom_names)
		atom_3_index = findfirst(x -> startswith(x, atom_3), target_atom_names)
		atom_1_coords = target_atom_coords[atom_1_index]
		atom_2_coords = target_atom_coords[atom_2_index]
		atom_3_coords = target_atom_coords[atom_3_index]
		vector_2_1 = atom_2_coords .- atom_1_coords
		vector_3_1 = atom_3_coords .- atom_1_coords
		distance_vector = cross_product_tuple(vector_2_1, vector_3_1)
		dot_product = dot(Electric_Field, distance_vector)
		angle = rad2deg(acos(dot_product / (norm(Electric_Field) * norm(distance_vector))))

		projection = (dot_product / norm(distance_vector)^2) .* distance_vector

		if angle > 90
			magnitude_projection = -norm(projection)
		elseif angle <= 90
			magnitude_projection = norm(projection)
		end

		electric_field_results[time] = Electric_Field
		angle_results[time] = angle
		projection_results[time] = projection
		magnitude_projection_results[time] = magnitude_projection
	end

	open("electric_field.log", "w") do io
		write(io, "Time,Electric_Field,Magnitude,Angle,Projection,Magnitude_Projection\n")

		for (time, electric_field) in electric_field_results
			magnitude = norm(electric_field)

			angle = angle_results[time]

			projection = projection_results[time]

			magnitude_projection = magnitude_projection_results[time]

			write(io, "$time,$electric_field,$magnitude,$angle,$projection,$magnitude_projection\n")
		end
	end

	angle_plot = scatter(xlabel = "t", ylabel = L"$\mathrm{\alpha}$", ylim=(0, 180))
	scatter!(sort(collect(keys(angle_results))), collect(values(sort(angle_results))), lw = 2, label=nothing)
	savefig(angle_plot, "angle_vs_time.pdf")

	cos_angle_results = Dict{Float64, Float64}()
	for (time, angle) in angle_results
		cos_angle_results[time] = cosd(angle)
	end

	cos_angle_plot = scatter(xlabel = "t", ylabel = L"cos($\mathrm{\alpha}$)")
	scatter!(sort(collect(keys(cos_angle_results))), collect(values(sort(cos_angle_results))), lw = 2, label=nothing)
	savefig(cos_angle_plot, "cos_angle_vs_time.pdf")

	proj_plot = scatter(xlabel = "t", ylabel = L"$\mathrm{scalar~projection~of~}\vec{p}$")
	scatter!(sort(collect(keys(magnitude_projection_results))), collect(values(sort(magnitude_projection_results))), lw = 2, label=nothing)
	savefig(proj_plot, "proj_vs_time.pdf")

	return electric_field_results, angle_results
end


function rdf_misc_all(dr, max_distance, iX, iY, iZ, fX, fY, fZ)
	results = process_gromacs_files()
	centre_of_mass_dict_unk = centre_of_mass_unk(results)
	centre_of_mass_dict_sol = centre_of_mass_sol(results)

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

	open("rdf_misc_all.log", "w") do io
		write(io, "Generating RDF plots...\n")
		for time in sort(collect(keys(centre_of_mass_dict_unk)))

			rdf_histogram = Dict{Float64, Float64}()
			for r in 0:dr:max_distance
				rdf_histogram[r] = 0.0
			end

			com_list_unk = centre_of_mass_dict_unk[time]
			com_list_sol = centre_of_mass_dict_sol[time]

			write(io, "Time: $time\n")

			for (unk_molecule_name, unk_centre_of_mass) in com_list_unk

				for (sol_molecule_name, com) in com_list_sol

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
			total_molecules = length(com_list_sol)
			for r in keys(time_step_rdf)
				time_step_rdf[r] /= (1 * total_molecules)
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

	num_snapshots = length(keys(centre_of_mass_dict_unk))
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

function electric_field_at_point_all(molecule_type, atom_1, atom_2, atom_3, iX, iY, iZ, fX, fY, fZ)
	k = 1.43996454 # Coulomb's constant in (eV nm^2 {elemental charge}^(-2))
	Lx = fX - iX
	Ly = fY - iY
	Lz = fZ - iZ

	atomic_charges = get_atomic_charges()

	electric_field_results = Dict{Float64, Dict{String, Tuple{Float64, Float64, Float64}}}()
	angle_results = Dict{Float64, Dict{String, Float64}}()
	projection_results = Dict{Float64, Dict{String, Tuple{Float64, Float64, Float64}}}()
	magnitude_projection_results = Dict{Float64, Dict{String, Float64}}()
	translations = [(dx*Lx, dy*Ly, dz*Lz) for dx in -1:1 for dy in -1:1 for dz in -1:1]
	results = process_gromacs_files()
	centre_of_mass_dict_unk = centre_of_mass_unk(results)

	for time in sort(collect(keys(results)))
		com_list_unk = centre_of_mass_dict_unk[time]

		for (unk_molecule_name, unk_centre_of_mass) in com_list_unk

			if molecule_type == "UNK" && endswith(unk_molecule_name, "UNL")
				continue
			elseif molecule_type == "UNL" && endswith(unk_molecule_name, "UNK")
				continue
			end

			middle_point = middle_of_atoms(results, unk_molecule_name, atom_1, atom_2, time)

			Electric_Field = (0.0, 0.0, 0.0)
			molecule_names, atom_names, atom_coords, box_dimensions = results[time]

			target_atom_names = atom_names
			target_atom_coords = atom_coords
			target_molecule_names = molecule_names

			for (molecule, atom_name, (x, y, z)) in zip(target_molecule_names, target_atom_names, target_atom_coords)

				if molecule == unk_molecule_name continue end

				q = nothing

				for key in keys(atomic_charges)
					if startswith(atom_name, key)
						q = atomic_charges[key]
						break
					end
				end

				if q === nothing
					continue
				end
				
				for translation in translations
					translated_atom_coords = (x, y, z) .+ translation
					r = translated_atom_coords .- middle_point
					d = norm(r)
					E_i = (k * q) .* r ./ d^3
					Electric_Field = Electric_Field .+ E_i
				end
			end
			
			atom_1_index = findfirst(x -> startswith(x, atom_1), target_atom_names)
			atom_2_index = findfirst(x -> startswith(x, atom_2), target_atom_names)
			atom_3_index = findfirst(x -> startswith(x, atom_3), target_atom_names)
			atom_1_coords = target_atom_coords[atom_1_index]
			atom_2_coords = target_atom_coords[atom_2_index]
			atom_3_coords = target_atom_coords[atom_3_index]
			vector_2_1 = atom_2_coords .- atom_1_coords
			vector_3_1 = atom_3_coords .- atom_1_coords
			distance_vector = cross_product_tuple(vector_2_1, vector_3_1)
			dot_product = dot(Electric_Field, distance_vector)
			angle = rad2deg(acos(dot_product / (norm(Electric_Field) * norm(distance_vector))))

			projection = (dot_product / norm(distance_vector)^2) .* distance_vector

			if angle > 90
				magnitude_projection = -norm(projection)
			elseif angle <= 90
				magnitude_projection = norm(projection)
			end

			if haskey(electric_field_results, time)
				electric_field_results[time][unk_molecule_name] = Electric_Field
			else
				electric_field_results[time] = Dict(unk_molecule_name => Electric_Field)
			end

			if haskey(angle_results, time)
				angle_results[time][unk_molecule_name] = angle
			else
				angle_results[time] = Dict(unk_molecule_name => angle)
			end

			if haskey(projection_results, time)
				projection_results[time][unk_molecule_name] = projection
			else
				projection_results[time] = Dict(unk_molecule_name => projection)
			end

			if haskey(magnitude_projection_results, time)
				magnitude_projection_results[time][unk_molecule_name] = magnitude_projection
			else
				magnitude_projection_results[time] = Dict(unk_molecule_name => magnitude_projection)
			end
		end
	end
	return electric_field_results, angle_results, magnitude_projection_results
end

function count_angle_in_cuboid_z(centre_of_mass_unk_dict, angles_dict, iX, iY, iZ, fX, fY, fZ)
	count_angle_in_cuboid_z_dict = Dict{Float64, Vector{Int}}()

	for time in sort(collect(keys(centre_of_mass_unk_dict)))
		com_list = centre_of_mass_unk_dict[time]
		angle_count = zeros(Int, 181)
		com_indices = []

		for (idx, com) in com_list
 			x, y, z = com
 			if iX <= x <= fX && iY <= y <= fY && iZ <= z <= fZ
 				push!(com_indices, idx)
			end
		end

		angle_list = angles_dict[time]

		for molecule_name in com_indices
			if !haskey(angle_list, molecule_name)
				continue
			end
			angle = angle_list[molecule_name]
			angle_int = round(Int, angle)
			angle_count[angle_int + 1] += 1
		end

		count_angle_in_cuboid_z_dict[time] = angle_count
	end

	return count_angle_in_cuboid_z_dict
end

function collect_magnitude_projection_in_cuboid_z(centre_of_mass_unk_dict, magnitude_projection_dict, iX, iY, iZ, fX, fY, fZ)
	collect_magnitude_projection_in_cuboid_z_dict = Dict{Float64, Vector{Float64}}()

	for time in sort(collect(keys(centre_of_mass_unk_dict)))
		com_list = centre_of_mass_unk_dict[time]
		magnitude_projection_collection = Float64[]
		com_indices = []

		for (idx, com) in com_list
 			x, y, z = com
 			if iX <= x <= fX && iY <= y <= fY && iZ <= z <= fZ
 				push!(com_indices, idx)
			end
		end

		magnitude_projection_list = magnitude_projection_dict[time]

		for molecule_name in com_indices
			if !haskey(magnitude_projection_list, molecule_name)
				continue
			end
			magnitude_projection = magnitude_projection_list[molecule_name]
			push!(magnitude_projection_collection, magnitude_projection)
		end

		collect_magnitude_projection_in_cuboid_z_dict[time] = magnitude_projection_collection
	end

	return collect_magnitude_projection_in_cuboid_z_dict
end

function bin_magnitude_projections(magnitude_projections, bin_edges)
	bin_counts = zeros(length(bin_edges) - 1)
	for projection in magnitude_projections
		for i in 1:(length(bin_edges) - 1)
			if bin_edges[i] <= projection < bin_edges[i + 1]
				bin_counts[i] += 1
				break
			end
		end
	end

	total_counts = sum(bin_counts)
   bin_counts_normalized = bin_counts / total_counts
   return bin_counts, bin_counts_normalized
end

function write_angle_in_multiple_slices(molecule_type, atom_1, atom_2, atom_3, iX, iY, iZ, fX, fY, fZ, Nz)

	results = process_gromacs_files()
	centre_of_mass_unk_dict = centre_of_mass_unk(results)
	electric_field_results, angle_results, magnitude_projection_results = electric_field_at_point_all(molecule_type, 
	atom_1, atom_2, atom_3, iX, iY, iZ, fX, fY, fZ)

    bar_data_avg = Dict()

    open("count_angle_in_multiple_slices.log", "w") do io
        println(io, "Counting molecules in multiple slices with angle populations...\n")
        println(io, "Initial point of the slices (iX, iY, iZ): ($iX, $iY, $iZ)")
        println(io, "Final point of the slices (fX, fY, fZ): ($fX, $fY, $fZ)\n")

        deltaZ = (fZ - iZ) / Nz

		  bin_edges = -7.0:0.1:7.0
	     bin_counts_avg = Dict()
		  bin_counts_normalized_avg = Dict()

        for n in 1:Nz
		  		local_iZ = iZ + (n - 1) * deltaZ
            local_fZ = iZ + n * deltaZ
            
				count_angle_in_cuboid_z_dict = count_angle_in_cuboid_z(centre_of_mass_unk_dict, angle_results,
            iX, iY, local_iZ,
            fX, fY, local_fZ)

				angle_count_avg = zeros(Float64, 181)

				collect_magnitude_projection_in_cuboid_z_dict = collect_magnitude_projection_in_cuboid_z(centre_of_mass_unk_dict,
				magnitude_projection_results, iX, iY, local_iZ, fX, fY, local_fZ)

				p_magnitude = plot(title = "", xlabel = "Time (ns)", ylabel =L"Scalar Projection ($\mathrm{eV~e^{-1}~nm^{-1}}$)", legend = false, ylim=(-9, 9))
				avg_projections = Dict()

				bin_counts_avg_time = zeros(length(bin_edges) - 1)
				bin_counts_normalized_avg_time = zeros(length(bin_edges) - 1)
				all_magnitude_projections = Float64[]

				for (time, magnitude_projections) in sort(collect_magnitude_projection_in_cuboid_z_dict)
					scatter!(p_magnitude, fill(time, length(magnitude_projections)), magnitude_projections, markeralpha=0.15, color=:blue)
					avg_projections[time] = mean(magnitude_projections)

					bin_counts, bin_counts_normalized = bin_magnitude_projections(magnitude_projections, bin_edges)
	            bin_counts_avg_time .+= bin_counts
					bin_counts_normalized_avg_time .+= bin_counts_normalized

					append!(all_magnitude_projections, magnitude_projections)

					df = DataFrame(Time = time, Counts = bin_counts, Normalized_Counts = bin_counts_normalized)
	            CSV.write("bin_counts_slice_$(n)_time_$(time).csv", df)

				end

				if !isempty(avg_projections)
					df_avg_proj = DataFrame(Time = collect(keys(avg_projections)), Avg_Projection = collect(values(avg_projections)))
					CSV.write("avg_projections_slice_$(n).csv", df_avg_proj)
				end

				if !isempty(collect_magnitude_projection_in_cuboid_z_dict)
					all_times = Int[]
					all_magnitude_projections = Float64[]

					for (time, magnitude_projections) in sort(collect_magnitude_projection_in_cuboid_z_dict)
						append!(all_times, fill(time, length(magnitude_projections)))
						append!(all_magnitude_projections, magnitude_projections)
					end

    				df_all_mag_proj = DataFrame(Time = all_times, Magnitude_Projections = all_magnitude_projections)
    				CSV.write("all_magnitude_projections_slice_$(n).csv", df_all_mag_proj)
				end

				if !isempty(all_magnitude_projections)
    				kde_estimate = kde(all_magnitude_projections)

    				p_kde = plot(kde_estimate.x, kde_estimate.density,
    				title = "Smooth KDE of Magnitude Projections for Slice: $n",
    				xlabel = "Magnitude Projection", ylabel = "Density", linewidth=2, color=:blue, legend=false)

    				savefig(p_kde, "smooth_kde_slice_$(n).pdf")

					df_kde = DataFrame(Values = kde_estimate.x, Density = kde_estimate.density)
	            CSV.write("kde_slice_$(n).csv", df_kde)
				else
    				println("Skipping KDE calculation for slice $n as there are no magnitude projections.")
				end

	        	bin_counts_avg_time ./= length(collect_magnitude_projection_in_cuboid_z_dict)
				bin_counts_normalized_avg_time ./= length(collect_magnitude_projection_in_cuboid_z_dict)
	        	bin_counts_avg[n] = bin_counts_avg_time
				bin_counts_normalized_avg[n] = bin_counts_normalized_avg_time

	         p_histogram = bar(bin_edges[1:end-1], bin_counts_avg_time, 
            title = "Histogram of Magnitude Projections for Slice: $n",
            xlabel = "Magnitude Projection", ylabel = "Count", legend = false)

        		savefig(p_histogram, "magnitude_projection_histogram_slice_$(n).pdf")

				p_histogram_normalized = bar(bin_edges[1:end-1], bin_counts_normalized_avg_time,
		      title = "Normalized Histogram of Magnitude Projections for Slice: $n",
        		xlabel = "Magnitude Projection", ylabel = "Density", legend = false)

        		savefig(p_histogram_normalized, "normalized_magnitude_projection_histogram_slice_$(n).pdf")

				for (time, avg_projection) in sort(avg_projections)
					scatter!(p_magnitude, [time], [avg_projection], color=:red, markershape=:circle)
				end

            savefig(p_magnitude, "magnitude_projection_slice_$(n).pdf")

            println(io, "Slice $n\n")
            println(io, "Initial point of the slice (iX, iY, iZ): ($iX, $iY, $local_iZ)")
            println(io, "Final point of the slice (fX, fY, fZ): ($fX, $fY, $local_fZ)\n")

            for time in sort(collect(keys(count_angle_in_cuboid_z_dict)))
                angle_counts = count_angle_in_cuboid_z_dict[time]
					 angle_count_avg .+= angle_counts
                println(io, "Time: $time, angle counts in slice $n: $angle_counts")
            end

				angle_count_avg ./= length(count_angle_in_cuboid_z_dict)
				bar_data_avg[n] = angle_count_avg

        		total_count = sum(values(count_angle_in_cuboid_z_dict))
            num_time_steps = length(count_angle_in_cuboid_z_dict)
            average_count = total_count / num_time_steps
            println(io, "\nAverage number of molecules in slice $n over time: $average_count\n")
        end

        println(io, "End of files, counting complete.")
    end

	 for (slice, angle_count_avg) in bar_data_avg
	 	p_avg = bar(0:180, angle_count_avg, title = "Slice: $slice",
		xlabel = L"angle", ylabel = "Average Count", legend = false)
		savefig(p_avg, "angle_slice_avg_$(slice).pdf")
	 end
end

