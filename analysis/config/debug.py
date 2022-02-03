import json

with open("selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

product = jsonObject['GOODMUON']
overall = jsonObject['triggers']

print(product)
print('------')
print('------')

for trigger in overall:
    if trigger['name'] =='isVBF' and trigger['year'] ==2018: print(trigger['definition'])
