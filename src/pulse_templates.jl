abstract type PulseShape{T} end
struct GumbelPulse{T} <: PulseShape{T}
    sigma::T
    amplitude::T
end

make_pulse_dist(p::T) where {T <:PulseShape} = error("not implemented")
make_pulse_dist(p::GumbelPulse{T}) where {T} = mu -> Gumbel(mu, p.sigma)
