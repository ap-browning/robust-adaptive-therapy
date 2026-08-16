#################################
## DOWNLOAD DATA
#################################

# Download data if it hasn't been downloaded already
if length(readdir("$(@__DIR__)/data")) < 100
    @warn("Downloading data from GitHub...")
    url = "https://github.com/reneebrady/IADT_PCaSC/raw/refs/heads/master/Archive.zip"
    tmp = download(url)
    dir = "$(@__DIR__)/data"
    zipreader = ZipFile.Reader(tmp)
    for f in zipreader.files
        if startswith(f.name,"1_SensitivityAnalysis") && endswith(f.name,"txt")
            # Full path
            path = joinpath(dir,split(f.name,"/")[end])
            # Create directory if needed
            if endswith(f.name,"/")
                mkdir(path)
            else
                write(path,read(f))
                # Processing... remove spaces, formatting issues, etc...
                content = read(path,String)
                content = replace(content," " => "")
                content = replace(content,"\""  => "")
                # Write file
                open(path,"w") do ff
                    write(ff, content)
                end
            end
        end
    end
    close(zipreader)
    rm(tmp)
end


#################################
## LOAD DATA FOR A PATIENT
#################################

# Patient list... (all)
fnames = filter(x -> split(x,".")[2] == "txt",readdir("$(@__DIR__)/data"))
patient_list = sort([parse(Int,split(f,".")[1][2:end]) for f in fnames])


"""
    Load data associated with a patient
"""
function patient_data(idx;normalised=true,output=:nice)
    idx ∈ patient_list || error("Patient $idx doesn't exist!")

    # Load
    data = CSV.read("$(@__DIR__)/data/p$idx.txt", DataFrame)

    # Output raw if requested
    if output == :raw
        return data
    end

    # Load treatment schedule
    datatmp = dropmissing(data[:,[:Day,:Treatment]])
    treat = datatmp.Treatment
    switch = [-Inf; datatmp.Day[findall(treat[1:end-1] .!= treat[2:end])]]
    Tx = t -> isodd(findlast(t .> switch))

    # Remove missing...
    data = dropmissing(data[:,[:Day, :PSA]])

    day = Int.(data.Day)
    psa = data.PSA
    if normalised
        psa = data.PSA / data.PSA[1]
    end

    # Remove missing
    filt_idx = psa .> 0.0
    day = Float64.(day[filt_idx])
    psa = Float64.(psa[filt_idx])

    return day, psa, Tx

end


# Patient list (minimum data)
patient_list_mindata_idx = falses(length(patient_list))
for i = eachindex(patient_list)
    day, psa, Tx = patient_data(patient_list[i])
    num_on = count(Tx.(day) .== true)
    num_off = length(day) - num_on
    patient_list_mindata_idx[i] = num_on ≥ 5 && num_off ≥ 5
end
patient_list_mindata = patient_list[patient_list_mindata_idx]