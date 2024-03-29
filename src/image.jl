"""
    image(x,sides)

Move coordinates to minimum periodic image.

"""
image(x, sides) = image_one.(x, sides)

# Translate a single coordinate according to periodic conditions
function image_one(x, side)
    xmod = mod(x, side)
    if xmod <= -side / 2
        return xmod + side
    elseif xmod > side / 2
        return xmod - side
    else
        return xmod
    end
end

