# find julia files to be converted into README.md

function gernerate_README(examplestoconvert, folder)

    ps = map(examplestoconvert) do str
        fs = filter(s -> split(s, ".")[end] == "jl", readdir(folder*str))
        if length(fs) > 0
            [str, fs[1]]
        else
            String[]
        end
    end
    filter!(p -> p != String[], ps)

    # This code below writes the README.md
    for p in ps
        s1 = string(folder,p[1],"/",p[2])
        str_arr = split(read(open(s1,"r"), String), "\n")

        # convert to md format
        str_arr = map(str_arr) do s
            if length(s) > 2 && s[1] == '#' && s[2] == ' '
              string(s[3:end],"\n")
            else
              string(s,"\n")
            end
        end

        s2 = string(folder,p[1],"/README.md")
        f = open(s2,"w")
        write(f, string(str_arr...))
        close(f)
    end
end
