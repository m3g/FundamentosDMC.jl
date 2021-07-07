"""

```
printxyz(time,x,pars,buffer)
```

Subroutine that writes the current coordinates to an output file, previously opened as `buffer`

"""
function printxyz(time,x,opt,buffer)
  n = length(x)
  println(buffer,n)
  println(buffer," time = ", time)
  for i in 1:n
    xt = image(x[i],opt.sides)
    println(buffer,@sprintf("He %12.5f %12.5f %12.5f", xt[1], xt[2], 0.))
  end
end