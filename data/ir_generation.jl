using WAV

long_ir_filename = "default_IR_48.wav"

long_file, sr, _, _ = wavread(long_ir_filename)


end_curve = 0.5 .* cos.(LinRange(0, pi, 100)) .+ 0.5

short_file = long_file[begin : 1024]
short_file[end-99 : end] .*= end_curve

open("base_ir.h", "w") do io 

    write(io, "global_const u32 ir_size = 1024;\n\n")
    write(io, "global_const float IR[ir_size] = {")
    
    for element in short_file
        write(io, "$(string(Float32(element)))f, ")
    end
    
    write(io, "};\n")
end 


println("done")
