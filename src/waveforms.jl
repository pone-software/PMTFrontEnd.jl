
struct PMTWaveform{T}
    photon_times::AbstractVector{T}
    photon_charges::AbstractVector{T}
    pulse_shape::PulseShape{T}
end


function evaluate_waveform(times, wf::PMTWaveform{T}) where T
    pulse_dist = make_pulse_dist(wf.pulse_shape)
    
    evaluated_wf = zeros(T, size(times))

    for (ptime, pcharge) in zip(wf.photon_times, wf.photon_charges)
        evaluated_wf += pdf.(pulse_dist(ptime), times) .* pcharge .* wf.pulse_shape.amplitude
    end

    evaluated_wf

end

function add_gaussian_white_noise(values, scale)
    values .+ randn(size(values)) * scale
end


struct DigitizedWaveform{T}
    timestamps::AbstractVector{T}
    values::AbstractVector{T}
end


function digitize_waveform(
    true_waveform::PMTWaveform{T},
    sampling_frequency::T,
    digitizer_frequency::T,
    noise_amp::T,
    filter,
    eval_range::Tuple{T, T},
    ) where T
    
    dt = 1/sampling_frequency # ns
    timesteps = range(eval_range[1], eval_range[2], step=dt)

    waveform_values = evaluate_waveform(timesteps, true_waveform)
    waveform_values_noise = add_gaussian_white_noise(waveform_values, noise_amp)

    waveform_filtered = filt(filter, waveform_values_noise)

    resampling_rate = digitizer_frequency / sampling_frequency
    new_interval = range(eval_range[1], eval_range[2], step=1/digitizer_frequency)
    waveform_resampled = resample(waveform_filtered, resampling_rate)

    new_interval, waveform_resampled
end

function make_nnls_matrix(
    pulse_times::AbstractVector{T},
    pulse_shape::PulseShape,
    timestamps::AbstractVector{T}) where {T}

    nnls_matrix = zeros(T, size(timestamps, 1), size(pulse_times, 1))

    for i in 1:size(pulse_times, 1)
        pulse = make_pulse_dist(pulse_shape)(pulse_times[i])
        nnls_matrix[:, i] = pulse_shape.amplitude .* pdf.(pulse, timestamps)
    end

    nnls_matrix
    
end


function apply_nnls(
    pulse_times::AbstractVector{T},
    pulse_shape::PulseShape{T},
    timestamps::AbstractVector{T},
    waveform::AbstractVector{T}) where {T}

    matrix = make_nnls_matrix(pulse_times, pulse_shape, timestamps)
    charges = nonneg_lsq(matrix, waveform; alg=:nnls)[:, 1]
end


