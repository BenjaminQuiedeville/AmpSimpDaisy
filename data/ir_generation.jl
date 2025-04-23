using WAV
using DSP
using Statistics

long_ir_filename = "default_IR_48.wav"
# long_ir_filename = "D:/Projets musique/ImpulseResponses/Rainbows/48/02 DV30.wav"


long_file, sr, _, _ = wavread(long_ir_filename)

short_size = 1024

tail_env_size = 30
tail_env = 0.5 .* cos.(LinRange(0, pi, tail_env_size)) .+ 0.5

short_file = long_file[begin : short_size]
short_file[end-(tail_env_size-1) : end] .*= tail_env

short_file .-= mean(short_file)

wavwrite(short_file, "base_ir.wav"; Fs = sr)

# short_file = zeros(Float64, short_size);
# short_file[1] = 1.0

# open("base_ir.h", "w") do io 

#     block_size = 8
#     part_size = 4 * block_size
#     nparts  =  Int(length(short_file) / block_size)

#     reshaped = reshape(short_file, block_size, :)
    
#     write(io, "global_const float ir_parts[npartitions][fft_size] = {\n")
    
#     for i in 1:nparts
        
#         time_signal = zeros(part_size)
#         time_signal[1:block_size] .= reshaped[:, i]
        
#         dft = DSP.rfft(time_signal)
        
#         write(io, "{\n\t")
        
#         for bin_index in 1:length(dft)
#             write(io, "$(string(Float32(dft[bin_index].re)))f, ")
#         end
#         write(io, "\n\t")
        
#         for bin_index in 2:(length(dft)-1)
#             write(io, "$(string(Float32(dft[bin_index].im)))f, ")
#         end
    
#         write(io, "\n},\n")
#     end
    
#     write(io, "};\n")
# end

# open("base_ir.h", "w") do io 

    # dft = DSP.rfft(short_file)
    # write(io, "alignas(16) global_const float ir_dft[dft_buffer_size] = {\n")
    
    # for index in 1:Int(fft_size/2)
    #     write(io, "$(string(Float32(dft[index].re)))f, ")
    # end
    
    # write(io, "\n")
    
    # for index in 1:Int(fft_size/2) - 1
    #     write(io, "$(string(Float32(dft[index].im)))f, ")
    # end
# end 

open("base_ir.h", "w") do io 
    
    write(io, "global_const u16 ir_size = $(short_size);\n")

    write(io, "global_const float ir[ir_size] = {\n\t")
    
    for elem in short_file
        write(io, "$(string(Float32(elem)))f, ")
    end
    
    write(io, "\n};\n")
end 

println("done")
