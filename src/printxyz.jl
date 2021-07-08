"""

```
printxyz(time,x,sys,buffer)
```

Subroutine that writes the current coordinates to an output file, previously opened as `buffer`.

"""
function printxyz(time,x,sys,buffer)
  n = length(x)
  println(buffer,n)
  println(buffer," time = ", time)
  for i in 1:n
    xt = image(x[i],sys.sides)
    if length(xt) == 2
      println(buffer,@sprintf("He %12.5f %12.5f %12.5f", xt[1], xt[2], 0.))
    else
      println(buffer,@sprintf("He %12.5f %12.5f %12.5f", xt[1], xt[2], xt[3]))
    end
  end
end

"""

```
printxyz(x,sys,file::String) 
```

Prints one set of coordinates, given in `x` to file `file`. Overwrites
the file.  

"""
function printxyz(x,sys,file::String) 
  buffer = open(file,"w") 
  printxyz(0,x,sys,buffer)
  close(buffer)
  return nothing
end