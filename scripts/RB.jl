using Cliffords, LsqFit

function clean(m::Void)
    return ""
end

function clean(m::SubString)
    return m
end

"""
Parses a GST experiment sequence. The first argument is a string describing the
    experiment sequence, while the second argument is a dictionary mapping gate
    labels to integers corresponding to indices in a pulse library. The output is
    an array of integers corresponding to indices in a pulse library.

    Experiment sequences have the following format. Gates are denoted by string 
    matching the regular expression "G[^G]+", e.g., 

       "G1", "Ga", or "G\$"

    A gate sequence matches a sequence of gates as defined above,  e.g., 

       "G1GbG#GaG%Gg" 

    An experiment sequence matches a sequence of gates, potentially followed by 
    a sequence power, i.e., a gate sequence in parenthesis followed by a carret
    and a non-negative integer (indicating repetition), followed
    by at most another gate sequence, e.g.,
  
       "(GaG1G\$)^4"

       "G1Gb(GcG]\$)^4G5G10"

       "G1Gb(GcG\$)^4"

       "(GcG\$)^4G5G10"

       "G1GbG5G10"
    

"""
function parse_gst(str, dict)
    function parse_count(str)
        count = r"\^([0-9]*).*"
        if length(str)>1 && str[1] == '^'
            parse(Int,match(count,str).captures[1])
        else
            1
        end
    end

    if contains("{}",str)
        return [dict["{}"]]
    else
        # split on parenthesis
        parts = split(str,['(',')'])
        # check if there is repetition (in the form of "^N")
        counts = map(parse_count,parts)
        counts = vcat(counts[2:end],[counts[1]])
        # extract subsequence strings
        clean_parts = [ match(r"(\^[0-9]*|)((G[^G]+)*)",part).captures[2] for part in parts ] 
        # map subsequences to arrays of gate labels
        clean_arr = map(s->split(s,'G'),clean_parts)
        # discard empty sequences
        nontrivial_arr = Vector{Int}[]
        nontrivial_counts = Int[]
        for (index,item) in enumerate(clean_arr)
            if length(item) > 1
                push!(nontrivial_arr, map(s->dict[s],item[2:end]))
                push!(nontrivial_counts, counts[index])
            end
        end
        return vcat([repmat(nontrivial_arr[i], nontrivial_counts[i]) for i in 1:length(nontrivial_arr)]...)
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

"""
FitResult type that stores the output of a fit.
"""
struct FitResult
    fit_params::Dict{String,Float64}
    sq_error::Float64
    Nσ::Float64
    errors::Dict{String,Float64}
    fit_curve::Function
    model_str::String
end

"""
`generic_fit(xpts, ypts, model, initial_guess, fit_params, model_string; yvars=[])`
Helper function to return a fit to a model.
# Arguments
* xpts: Independent fit variable.
* ypts: Dependent fit variable.
* model: Fit function(x, p) of data x and parameter vector p.
* initial_guess: Initial guess of fit parameters.
* fit_params::function(Array{Float64} => Dict{String, Float64}): Result dict of fit parameters.
* model_string::String: String indicating fit function.
* yvars: Variance of y points. Defaults to empty vector.
# Returns
* A FitResult object for the specified fit.
"""
function generic_fit(xpts, ypts, model, initial_guess, fit_params, model_string::String; yvars=[])
    @assert length(xpts) == length(ypts) "X and Y data length must match."
    #try
        if isempty(yvars)
          result = curve_fit(model, xpts, ypts, initial_guess)
          errors = estimate_errors(result)
          sq_error = sum(((model(xpts, result.param) - ypts) .^ 2))
        else
          @assert length(ypts) == length(yvars) "Y data and Y variance lengths must match."
          result = curve_fit(model, xpts, ypts, 1 ./ sqrt(yvars), initial_guess)
          errors = estimate_errors(result)
          sq_error = sum(((model(xpts, result.param) - ypts) .^ 2) ./ yvars)
        end
    #catch
    #    result =
    # Compute badness of fit:
    # Under the null hypothesis (that the model is valid and that the observations
    # do indeed have Gaussian statistics), the mean squared error is χ² distributed
    # with `dof` degrees of freedom. We can quantify badness-of-fit in terms of how
    # far the observed MSE is from the expected value, in units of σ = 2dof (the expected
    # standard deviation for the χ² distribution)
    dof = length(xpts)-length(initial_guess)
    Nσ = sq_error/sqrt(2*dof) - dof/sqrt(2*dof)
    return FitResult( fit_params(result.param),
              sq_error,
              Nσ,
              fit_params(errors),
              xpts->model(xpts,result.param),
              model_string)
end

"""`fit_RB(xpts, ypts, yvars=[])`
Fit to RB decay a*(1-b)^x + c
"""
function fit_RB(xpts, ypts; yvars=[])
    RB_fit_dict(p) = Dict("a" => p[1], "b" => p[2], "c" => p[3])
    model(n, p) = p[1] * (1-p[2]).^n + p[3]
    p_guess = [0.5, .01, 0.5]
    return generic_fit(xpts, ypts, model, p_guess, RB_fit_dict, "a * exp(-n/b) + c", yvars=yvars)
end
