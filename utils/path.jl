#### Paths ####

struct Path
    src
    dst
    edges
    cost
end

function Path()
    return Path(1, 1, [], 10000)
end

