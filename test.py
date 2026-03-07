import csv
with open('trialData.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Trial', 'Angular Diameter', 'Alpha Value'])

with open('trialData.csv', 'a', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['4', '1.52', '0.13'])
    writer.writerow(['5', '1.5', '.1'])