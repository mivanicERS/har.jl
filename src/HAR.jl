module HAR

using NamedArrays

readHar = function (fileName, useCoefficientsAsNames = false,
           toLowerCase = true)

    # HAR files come from the age of paper tapes--we need to get records first
    records = getRecords(fileName)

    # Some records define the beginnings of headers
    headers = Dict{String,Array}()
    currentHeader = ""

    for rr = records
        if length(rr) == 4
            currentHeader = read(IOBuffer(rr), String)

            # Headers may have trailing spaces on the right; remove them
            currentHeader = rstrip(currentHeader)

            # If we are doing lower case, make sure header names are lower case too
            if toLowerCase
                currentHeader = lowercase(currentHeader)
            end

            headers[currentHeader] = []
        else
            push!(headers[currentHeader], rr)
        end
    end

    # We return a dictionary with each header's data
    toRet = Dict()

    for hh in headers
        headerContent = processHeader(hh[2], toLowerCase)
        nameHeader = useCoefficientsAsNames ? (haskey(toRet,"coefficient") ? toRet["coefficient"] : hh[1]) : hh[1]
        toRet[nameHeader] = headerContent["values"]
    end



    return toRet
end

getRecords = function (fileName)

    bf = Base.open(fileName)

    records = []

    while !Base.eof(bf)

        recordSize = read(bf, Int32)

        recordData = Vector{UInt8}()
        for rs = 1:recordSize
            push!(recordData, read(bf, UInt8))
        end

        recordSizeCheck = read(bf, Int32)

        if recordSizeCheck != recordSize
            printl("Inconsistent record size")
        end

        push!(records, recordData)

    end


    Base.close(bf)

    return records
end

processHeader = function (hh, toLowerCase)

    toRet = Dict()

    toRet["type"] = read(IOBuffer(hh[1][5:10]), String)
    toRet["description"] = read(IOBuffer(hh[1][11:80]), String)
    toRet["numberOfDimensions"] = read(IOBuffer(hh[1][81:84]), Int32)

    toRet["dimensions"] = map(i -> read(IOBuffer(hh[1][(85+(i-1)*4):(84+i*4)]), Int32), 1:toRet["numberOfDimensions"])

    if toRet["type"] == "1CFULL"

        combinedData = reduce((a, f) -> append!(f[17:length(f)], a), hh[2:length(hh)], init=Vector{UInt8}())

        combinedString = read(IOBuffer(combinedData), String)


        toRet["values"] = map(f -> combinedString[((f-1)*toRet["dimensions"][2]+1):(f*toRet["dimensions"][2])] |> strip, 1:toRet["dimensions"][1])

        if toLowerCase == true
            toRet["values"] = map(lowercase, toRet["values"])
        end

    elseif toRet["type"] == "2IFULL"
        combinedData = reduce((a, f) -> append!(f[17:length(f)], a), hh[2:length(hh)], init=Vector{UInt8}())
        combinedBuffer = IOBuffer(combinedData)
        combinedIntegers = map(f -> read(combinedBuffer, Int32), 1:(toRet["dimensions"][1]*toRet["dimensions"][2]))
        toRet["values"] = reshape(combinedIntegers, (toRet["dimensions"][1], toRet["dimensions"][2]))
    elseif in(toRet["type"], ["REFULL", "RESPSE"])
        definedDimensions = read(IOBuffer(hh[2][5:8]), Int32)
        usedDimensions = read(IOBuffer(hh[2][13:16]), Int32)
        coefficient = strip(read(IOBuffer(hh[2][17:28]), String))

        toRet["coefficient"] = coefficient

        # Default dimensions
        allDimensions = ["dimensions"]
        dnames = 1:toRet["numberOfDimensions"]
        dimNames = Dict()
        for d = 1:toRet["numberOfDimensions"]
            dimNames[d] = 1:1:toRet["dimensions"][d]
        end
        uniqueDimNames = []
        if usedDimensions > 0
            allDimensions = read(IOBuffer(hh[2][33:(33+usedDimensions*12-1)]), String)
            dnames = map(f -> strip(allDimensions[((f-1)*12+1):(f*12)]), 1:usedDimensions)

            if toLowerCase
                dnames = map(lowercase,dnames)
            end

            dimNames = Dict()

            uniqueDimNames = unique(dnames)

            for d = 1:length(uniqueDimNames)
                nele = read(IOBuffer(hh[2+d][13:16]), Int32)
                allDim = read(IOBuffer(hh[2+d][17:(17+nele*12-1)]), String)
                dimNames[uniqueDimNames[d]] = map(f -> String(strip(allDim[((f-1)*12+1):(f*12)])), 1:nele)

                if toLowerCase == true
                    dimNames[uniqueDimNames[d]] = map(lowercase, dimNames[uniqueDimNames[d]])
                end

            end

        end

        dataStart = 2 + length(uniqueDimNames) + 1

        if toRet["type"] == "REFULL"
            numberOfFrames = read(IOBuffer(hh[dataStart][5:8]), Int32)

            numberOfDataFrames = convert(Int32, (numberOfFrames - 1) / 2)


            dataFrames = (dataStart) .+ (1:numberOfDataFrames) .* 2

            dataBytes = reduce((a, f) -> append!(a, hh[f][9:length(hh[f])]), dataFrames, init=Vector{UInt8}())

            numberOfValues = 1
            for (key, value) in dimNames
                numberOfValues = numberOfValues * length(value)
            end

            dims = map(f -> length(dimNames[f]), dnames)
            namesElements = map(f -> dimNames[f], dnames)

            numberOfValues = reduce((a, f) -> a * length(dimNames[f]), dnames, init=1)

            dataBytesBuffer = IOBuffer(dataBytes)

            dataVector = map(f -> read(dataBytesBuffer, Float32), 1:numberOfValues)


            toRet["values"] = NamedArray(reshape(dataVector, Tuple(dims)), Tuple(namesElements), Tuple(dnames))

        else
            elements = read(IOBuffer(hh[dataStart][5:8]), Int32)
            numberOfValues = reduce((a, f) -> a * length(dimNames[f]), dnames, init=1)
            dataVector = zeros(numberOfValues)

            for rr = (dataStart+1):length(hh)

                dataBytes = hh[rr][17:length(hh[rr])]

                currentPoints = convert(Int64, length(dataBytes) / 8)


                locationsBuffer = IOBuffer(dataBytes[1:(4*currentPoints)])

                locations = map(f -> read(locationsBuffer, Int32), 1:currentPoints) #readBin(dataBytes[1:(4*currentPoints)],'integer',size=4, n = currentPoints)

                valuesBuffer = IOBuffer(dataBytes[(4*currentPoints+1):(8*currentPoints)])

                values = map(f -> read(valuesBuffer, Float32), 1:currentPoints) #readBin(dataBytes[(4*currentPoints+1):(8*currentPoints)],'double',size=4, n = currentPoints)

                dataVector[locations] = values

                dims = map(f -> length(dimNames[f]), dnames)
                namesElements = map(f -> dimNames[f], dnames)

                toRet["values"] = NamedArray(reshape(dataVector, Tuple(dims)), Tuple(namesElements), Tuple(dnames))
            end

        end
        # else

        #     # m = array(readBin(
        #     #     headers[[h]]$records[[length(headers[[h]]$records)]][9:length(headers[[h]]$records[[3]])],
        #     #     'double',
        #     #     size = 4,
        #     #     n = prod(headers[[h]]$dimensions)
        #     #   ),
        #     #   dim = headers[[h]]$dimensions)
        #     toRet["values"] = nothing
        # end
    else
        toRet["values"] = nothing
    end
    return toRet
end
end
