"""

```
printxyz(time,x,pars,buffer)
```

Subroutine that writes the current coordinates to an output file, previously opened as `buffer`

"""
function printxyz(time,x,pars,buffer)
  n = length(x)
  println(buffer,n)
  println(buffer," time = ", time)
  for i in 1:n
    xt = image(x,pars.sides)
    println(buffer,"H $(xt[1]) $(xt[2]) 0.")
  end
end