using Cliffords, Formatting

"""
Clifford circuit randomization and compilation with recovery operation.

Since the recovery operation is included, there is no mask to output
(it essentially corresponds to the recovery operation.)
"""
function compare(pseqs,rseqs)

    maxlen = maximum(map(length,pseqs))
    rmaxlen = maximum(map(length,rseqs))

    seq   = Vector{Clifford}(maxlen)
    rseq   = Vector{Clifford}(rmaxlen)

    # now on to experiments with non-zero length (that is why we start at index 2)
    for n in 1:length(pseqs)
        # accumulate sequences
        inter = 1 # identity index
        for i in 1:length(pseqs[n])
            inter   = multiply(pseqs[n][i],inter)
        end
        rinter = 1
        for i in 1:length(rseqs[n])
            rinter  = multiply(rseqs[n][i],rinter)
        end
        # check that sequences are same up to final random Z90 rot
        index = multiply(inter,Cliffords.localcliffordindex(inv(localclifford(rinter))))
        if !(index in [1,9])
            error("Mismatch in line $(n)")
        else
            index == 1 ? (@printf "I") : (@printf "Z")
        end

    end
end

begin
    global multiply

    table = Dict{Tuple{Int,Int},Int}()

    function multiply(c1,c2)
        c1_ = round(Int,c1==0 ? 1 : c1)
        c2_ = round(Int,c1==0 ? 1 : c2)
        if !((c1_,c2_) in keys(table))
            table[(c1_,c2_)] = Cliffords.localcliffordindex(localclifford(c1_)*localclifford(c2_))
        end
        return table[(c1_,c2_)]
    end
end

# begin
#     println("Comparing gst-1q-1024.csv and $(ARGS[1])")
#
#     seqs = readlines(open("gst-1q-1024.csv"))
#     pseqs = [map(x->parse(Int,x), split(s,',')) for s in seqs];
#     seqs = readlines(open(ARGS[1]))
#     rseqs = [map(x->parse(Int,x), split(s,',')) for s in seqs];
#
#     compare(pseqs,rseqs)
#     @printf "\n"
#     println("Done!")
# end

begin
    infile = open(ARGS[1])
    seqs = readlines(infile)
    close(infile)

    #range = parse(Int,ARGS[2])
    #count = parse(Int,ARGS[3])

    pseqs = [map(x->parse(Int,x), split(s,',')) for s in seqs];

    #L = ceil(Int,log(10,range)+1)
    for k in 1:10
        # i = rand(1:range)
        label = format("{1:02d}",k)
        println("Comparing $(ARGS[1]) and $(replace(ARGS[1],r".csv$" => "-r$(label).csv"))")
        rfile = open(replace(ARGS[1],r".csv$" => "-r$(label).csv"))
        rseqs = [map(x->parse(Int,x), split(s,',')) for s in readlines(rfile)];
        compare(pseqs, rseqs)
        @printf "\n"
        close(rfile)
    end
end
