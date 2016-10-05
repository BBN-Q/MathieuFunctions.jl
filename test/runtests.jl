using Base.Test, MathieuFunctions

@test maximum(CharacteristicA(0,k=0:100) - [0:100;].^2) == 0

@test norm(CharacteristicB(0,k=1:100) - [1:100;].^2) == 0

begin
    test1 = readcsv("./test/MathieuCharacteristicA-1.csv");
    @test (test1 - reduce(hcat,Vector{Float64}[CharacteristicA(q,k=0:10) for q in [-10:.01:10;]]) |> abs |> maximum) < 7.5e-13
end

begin
    test1 = readcsv("./test/MathieuCharacteristicA-2.csv");
    @test (test1 - reduce(hcat,Vector{Float64}[CharacteristicA(q,k=0:3) for q in [30:.01:50;]]) |> abs |> maximum) < 7.6e-13
end

begin
    test1 = readcsv("./test/MathieuCharacteristicB-1.csv");
    @test (test1 - reduce(hcat,Vector{Float64}[CharacteristicB(q,k=1:10) for q in [-10:.01:10;]]) |> abs |> maximum) < 7.5e-13
end

begin
    test1 = readcsv("./test/MathieuCharacteristicB-2.csv");
    @test (test1 - reduce(hcat,Vector{Float64}[CharacteristicB(q,k=1:3) for q in [30:.01:50;]]) |> abs |> maximum) < 2.8e-11
end

