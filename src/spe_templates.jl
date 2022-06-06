# Structs
abstract type SPEDistribution{T} end
struct ExponTruncNormalSPE{T}
    expon_rate::T
    norm_sigma::T
    norm_mu::T
    trunc_low::T
    expon_weight::T
end

# Distribution factories
make_spe_dist(d::SPEDistribution{T}) where T = error("not implemented")

function make_spe_dist(d::ExponTruncNormalSPE{T}) where T
    norm = Normal(d.norm_mu, d.norm_sigma)
    tnorm = truncated(norm, lower=d.trunc_low)

    expon = Exponential(d.expon_rate)
    dist = MixtureModel([expon, tnorm], [d.expon_weight, 1 - d.expon_weight])
        
    return dist

end