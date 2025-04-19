using WAV
using DSP

long_ir_filename = "default_IR_48.wav"

long_file, sr, _, _ = wavread(long_ir_filename)

short_ir_size = 512
dft_n_bins = 1024 / 2 + 1

end_curve = 0.5 .* cos.(LinRange(0, pi, 100)) .+ 0.5

# short_file = long_file[begin : 1024]
# short_file[end-99 : end] .*= end_curve

short_file = zeros(Float64, short_ir_size);
short_file[1] = 1.0

dft = DSP.rfft(short_file)


open("base_ir.h", "w") do io 

    write(io, "alignas(16) global_const float ir_dft[dft_buffer_size] = {\n")
    
    for element in dft
        write(io, "$(string(Float32(element.re)))f, ")
    end
    
    write(io, "\n")
    
    for element in dft
        write(io, "$(string(Float32(element.im)))f, ")
    end

    write(io, "\n};\n")
end 


println("done")
