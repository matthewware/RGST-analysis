using Cliffords, Formatting

function pauli_rand_clifford_circuit(pseqs, ofile)
    # Pauli ops 
    ps = [Id,X,Y,Z]
    # Pauli ops as Clifford group ops
    c_ps = map(localclifford,[1,3,6,9])

    for n in 1:length(pseqs)
        prop = Id # Pauli correction propagated through sequence
        rp   = Id
        for i in 1:length(pseqs[n])-1
            c = pseqs[n][i] == 0 ? localclifford(1) : localclifford(pseqs[n][i])
            # pick a new random pauli
            rp = rand(ps)
            # combine new Pauli with propagated Pauli to get new propagated Pauli
            # compute new "compiled" pulse
            prop, comp = propagate_and_compile(c,rp,prop)
            @printf ofile "%d," comp
        end
        c = pseqs[n][end] == 0 ? localclifford(1) : localclifford(pseqs[n][end])
        rp = rand(ps)
        prop, comp = propagate_and_compile(c,rp,prop)
        @printf ofile "%d\n" Cliffords.localcliffordindex(Clifford(rand([Id,Z])*prop)*(c*Clifford(rp)))
        #println(Cliffords.localcliffordindex(Clifford(rand([Id,Z])*prop)*(c*Clifford(rp))))
    end
end

function pauli_rand_clifford_circuit2(pseqs, ofile)
    # Pauli ops 
    paulis = [Id,X,Y,Z]
    # Pauli ops as Clifford group ops
    cliffs = [1,3,6,9]
    c_ps = map(localclifford,[1,3,6,9])

    for n in 1:length(pseqs)
        prop = Id # Pauli correction propagated through sequence
        rp   = Id
        for i in 1:length(pseqs[n])-1
            c = pseqs[n][i] == 0 ? localclifford(1) : localclifford(pseqs[n][i])
            # pick a new random pauli
            rp = rand(1:4)
            # combine new Pauli with propagated Pauli to get new propagated Pauli
            # compute new "compiled" pulse
            prop, _ = propagate_and_compile(c,paulis[rp],prop)
            @printf ofile "%d,%d," cliffs[rp] pseqs[n][i]
        end
        c = pseqs[n][end] == 0 ? localclifford(1) : localclifford(pseqs[n][end])
        rp = rand(1:4)
        prop, _ = propagate_and_compile(c,paulis[rp],prop)
        @printf ofile "%d,%d,%d\n" cliffs[rp] pseqs[n][end] Cliffords.localcliffordindex(Clifford(rand([Id,Z])*prop))
        # println(cliffs[rp],",",pseqs[n][end],",",Cliffords.localcliffordindex(Clifford(rand([Id,Z])*prop)))
    end
end

# As is, the clifford multiplication is still expensive, so we need to memoize it (cache it)
# in order to speed up the computation.
let
    global propagate_and_compile

    table = Dict{Tuple{Clifford,Pauli,Pauli},Tuple{Pauli,Int}}()

    function propagate_and_compile(cliff,randp,propp)
        if !((cliff,randp,propp) in keys(table))
            prop = cliff*(randp*propp)
            table[(cliff,randp,propp)] = (prop,Cliffords.localcliffordindex(cliff*Clifford(randp)))
        end
        return table[(cliff,randp,propp)]
    end
end

begin 
    # set random seed (useful for testing)
    srand(0xcafe)

    infile = open(ARGS[1])
    seqs = readlines(infile)
    close(infile)

    count = parse(Int,ARGS[2])

    pseqs = [map(x->parse(Int,x), split(s,',')) for s in seqs];

    L = ceil(Int,log(10,count)+1)
    for i in 1:count
        label = format("{1:0$(L)d}",i)
        println(replace(ARGS[1],r".csv$","-r$(label).csv"))
        rfile = open(replace(ARGS[1],r".csv$","-r$(label).csv"),"w")
        pauli_rand_clifford_circuit(pseqs, rfile)
        close(rfile)
    end
end
