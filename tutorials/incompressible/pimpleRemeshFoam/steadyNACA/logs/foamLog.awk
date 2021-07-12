# Awk script for OpenFOAM log file extraction
BEGIN {
    Iteration=0
    resetCounters()
}

# Reset counters used for variable postfix
function resetCounters() {
    # Reset counters for 'Solving for ...'
    for (varName in subIter)
    {
        subIter[varName]=0
    }
}

# Extract value after columnSel
function extract(inLine,columnSel,outVar,a,b)
{
    a=index(inLine, columnSel)
    b=length(columnSel)
    split(substr(inLine, a+b),outVar)
    gsub("[,:]","",outVar[1])
}

# Iteration separator (increments 'Iteration')
/^[ \t]*Time = / {
    Iteration++
    resetCounters()
}

# Time extraction (sets 'Time')
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    Time=val[1]
}

# Skip whole line with singularity variable
/solution singularity/ {
    next;
}

# Extract: 'Solving for ...'
/Solving for/ {
    extract($0, "Solving for ", varNameVal)

    varName=varNameVal[1]
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "Initial residual = ", val)
    print Time "\t" val[1] > file

    varName=varNameVal[1] "FinalRes"
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "Final residual = ", val)
    print Time "\t" val[1] > file

    varName=varNameVal[1] "Iters"
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "No Iterations ", val)
    print Time "\t" val[1] > file
}

# End
