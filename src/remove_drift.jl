"""

```
remove_drift!(v)
```

Removes possible drift from velocities.

"""
function remove_drift!(v)
  vmean = mean(v)
  for i in 1:length(v)
    v[i] -= vmean
  end
end