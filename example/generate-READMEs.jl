# find julia files to be converted into README.md

examplestoconvert = ["random_particles", "random_particles_in_circle"]

pairs = map(examplestoconvert) do str
    fs = filter(s -> split(s, ".")[end] == "jl", readdir(str))
    if length(fs) > 0
        [str, fs[1]]
    else
        String[]
    end
end
filter!(p -> p != String[], pairs)

# This code below writes the README.md
for p in pairs
    s1 = string("../example/",p[1],"/",p[2])
    str_arr = split(read(open(s1,"r"), String), "\n")

    # convert to md format
    str_arr = map(str_arr) do s
        if length(s) > 2 && s[1] == '#' && s[2] == ' '
          string(s[3:end],"\n")
        else
          string(s,"\n")
        end
    end

    s2 = string("../example/",p[1],"/README.md")
    f = open(s2,"w")
    write(f, string(str_arr...))
    close(f)
end
