from itertools import groupby

# Input array
a =   [
    0.9981,
    0.992,
    0.9834,
    0.971,
    0.9682,
    0.953,
    0.953,
    0.945,
    0.919,
    0.919,
    0.919,
    0.9997,
    0.9947,
    1.0005,
    1.0006
    ]
# Process the array
result = []
temp = []

for key, group in groupby(a):
    count = len(list(group))
    if count > 1:  # Group has more than one element
        if temp:  # Add combined single elements before adding a new group
            result.append(f"[{', '.join(map(str, temp))}]")
            temp = []
        result.append(f"[{key}]*{count}")
    else:  # Group has a single element
        temp.append(key)

# Add any remaining single elements
if temp:
    result.append(f"[{', '.join(map(str, temp))}]")

# Join the results into the desired format
output = " + ".join(result)

# Print the output
print(output)
