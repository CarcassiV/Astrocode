import time
 
start_time = time.time()
for i in range(1000000):
    print(i)
end_time = time.time()
print(f"Time taken with print: {end_time - start_time} seconds")
 
start_time = time.time()
for i in range(1000000):
    pass
end_time = time.time()
print(f"Time taken without print: {end_time - start_time} seconds")