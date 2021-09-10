N = 100
START = Array{Float64}(undef,N,N)

function next(mesh, step)
    newmesh = copy(mesh)
    sx=size(mesh)[1]
    sy=size(mesh)[2]
    for (i,j) in zip(1:size1, 1:size2)
        newmesh[i,j] = mesh[i,j] + step*(
        mesh[(i+1)%sx,j] + mesh[(i-1)%,j] + mesh[i,j] + mesh[i,j] -
        4mesh[i,j])
    end
end
