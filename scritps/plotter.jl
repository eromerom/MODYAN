#!/usr/bin/env julia

using CSV, Plots, DataFrames, LaTeXStrings

file_labels = []
file_data = []

for i in 1:2:length(ARGS)
    file_name = ARGS[i]
    file_label = ARGS[i+1]
    push!(file_labels, file_label)
    data = CSV.read(file_name, DataFrame)
    push!(file_data, data)
end

plot()

for i in 1:length(file_data)
    data = file_data[i]
    file_label = file_labels[i]
    
    values = data[:, "Values"]
    density = data[:, "Density"]
    
    plot!(values, density, label=file_label)
end

xlabel!(L"Scalar Projection (eV e$^{-1}$ nm$^{-1}$)")
ylabel!("Probability Density")
title!("")

savefig("output.pdf")

display(plot)

