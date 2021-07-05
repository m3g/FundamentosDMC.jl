"""

```
image(x,sides)
```

Move coordinates to minimum perioddic images

"""
function image(x,sides)
  image = mod.(x,sides)
  if image < -sides/2
    image = image + sides
  elseif image > sides/2
    image = image - sides
  end
  return image
end