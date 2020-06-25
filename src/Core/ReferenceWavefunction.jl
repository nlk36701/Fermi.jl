abstract type AbstractReferenceWavefunction end
begin 
    mutable struct ReferenceWavefunction{T} <: AbstractReferenceWavefunction where T <: AbstractFloat
        refEnergy::T
        vnuc::T
        nocca::Int
        noccb::Int
        nvira::Int
        nvirb::Int
        basis::B where B <: Lints.BasisSet
        molecule::M where M <: Lints.Molecule
        Ca::Array{T,2} #AO->MO coefficients a
        Cb::Array{T,2} #AO->MO coefficients b
        S::Array{T,2}
        T::Array{T,2}
        V::Array{T,2}
        epsa::Array{T,1} #orbital eigenvalues a
        epsb::Array{T,1} #orbital eigenvalues b
        ERI::I where I <: Fermi.AbstractTensor#AO basis electron repulsion integrals
    end
end

function ReferenceWavefunction(basis,molecule,nocca,noccb)
    ReferenceWavefunction(basis,molecule,nocca,noccb,
                          Fermi.ComputeEnvironment.interconnect,
                          Fermi.ComputeEnvironment.communicator,
                          Fermi.ComputeEnvironment.accelerator)
end

function ReferenceWavefunction(basis,molecule,nocca,noccb,
                               interconnect::Fermi.Environments.No_IC,
                               communicator::Fermi.Environments.NoCommunicator,
                               accelerator::Fermi.Environments.NoAccelerator)
    sz = Lints.getsize(S_engine,basis)
    S_engine = Lints.OverlapEngine(nprim,l)
    T_engine = Lints.KineticEngine(nprim,l)
    V_engine = Lints.NuclearEngine(nprim,l,molecule)
    I_engines = []
    for i in 1:Threads.nthreads()
        push!(I_engines,Lints.ERIEngine(nprim,l))
    end
    S = zeros(sz,sz)
    Ca = zeros(size(S))
    Cb = zeros(size(S))
    T = zeros(sz,sz)
    V = zeros(sz,sz)
    I = zeros(sz,sz,sz,sz)
    Lints.make_2D(S,S_engine,basis)
    Lints.make_2D(T,T_engine,basis)
    Lints.make_2D(V,V_engine,basis)
    Lints.make_ERI(I,I_engines,basis)
    ReferenceWavefunction{Float64}(
                                   0.0,
                                   0.0,
                                   nocca,
                                   noccb,
                                   sz-nocca,
                                   sz-noccb,
                                   basis,
                                   molecule,
                                   Ca,
                                   Cb,
                                   S,
                                   T,
                                   V,
                                   zeros(Float64,sz),
                                   zeros(Float64,sz)
                                   I
                                  )

end
