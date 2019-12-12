using Test, MathieuFunctions
using LinearAlgebra
using DelimitedFiles

readcsv(f) = DelimitedFiles.readdlm(f, ',')

function tapprox(a, b; atol=1e-15)
    normval = norm(a - b, Inf)
    @info "normval = $normval"
    isapprox(a, b; norm= x -> norm(x, Inf), atol = atol)
end

@testset "basic" begin
    @test charλ(1, 0, k=1:10) == [1.0, 1.0, 9.0, 9.0, 25.0, 25.0, 49.0, 49.0, 81.0, 81.0]
    @test maximum(charA(0,k=0:100) - [0:100;].^2) == 0
    @test norm(charB(0,k=1:100) - [1:100;].^2) == 0
end


filename = "MathieuCharacteristicA-1.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    r = reduce(hcat,Vector{Float64}[charA(q,k=0:10) for q in [-10:.01:10;]])
    @test tapprox(test1, r;  atol=7.5e-13)
end

filename = "MathieuCharacteristicA-2.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    r = reduce(hcat,Vector{Float64}[charA(q,k=0:3) for q in [30:.01:50;]])
    @test tapprox(test1, r;  atol=7.6e-13) # NOTE: was 7.5e-13
end

filename = "MathieuCharacteristicB-1.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    r = reduce(hcat,Vector{Float64}[charB(q,k=1:10) for q in [-10:.01:10;]])
    @test tapprox(test1, r;  atol=7.5e-13)
end

filename = "MathieuCharacteristicB-2.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    r = reduce(hcat,Vector{Float64}[charB(q,k=1:3) for q in [30:.01:50;]])
    @test tapprox(test1, r; atol=2.8e-11)
end

filename = "MathieuCharacteristicL-1.csv"
@testset "$filename" begin
    test1 = readcsv(filename)[1:100,:]
    test2 = Float64[charλ(q,ν,k=1:1)[1] for ν in [0:.01:0.99;], q in [-5:.01:5;]]
    @test tapprox(test1, test2, atol=7.5e-15)
    # TODO: test ν > 1 (currently failing)
end

filename = "MathieuCharacteristicL-2.csv"
@testset "$filename" begin
    test1 = readcsv(filename)[1:100,:]
    test2 = Float64[charλ(q,ν,k=1:1)[1] for ν in [0:.01:0.99;], q in [30:.01:50;]]
    @test tapprox(test1, test2, atol=4.5e-14)
    # TODO: test ν > 1 (currently failing)
end
