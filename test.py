import json
dicc = {}
with open("heartbeat_search.json", "w") as f:
    json.dump(dicc, f)

with open("heartbeat_search.json", "rb") as f:
    dicd = json.load(f)

print(dicd)