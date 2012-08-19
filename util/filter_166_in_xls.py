import xlrd
import xlwt

from draw_roc import get_166_manual_groups

groups = get_166_manual_groups()
pdbs = []
for g in groups:
    for i in g:
        pdbs.append(i)

new_wb = xlwt.Workbook()
new_sh = new_wb.add_sheet("166 data")

wb = xlrd.open_workbook("antiginInfo-rich.xls")
sh = wb.sheet_by_index(0)
row_count = 0
for rownum in range(sh.nrows):
    row = sh.row_values(rownum)
    pdb = '_'.join(row[:2])
    if pdb in pdbs:
        for i in xrange(len(row)):
            new_sh.write(row_count,i,row[i])
        row_count += 1

print row_count
new_wb.save("antigen_info_166.xls")
