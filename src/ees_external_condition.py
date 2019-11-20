import csv
import datetime
import sys
sys.path.append('../')

from src.air import Air
from src.temperature import Temperature


class ExternalCondition:

    def __init__(self, path: str, region: str, reading_order, error_print):

        column_num = {
            'region1': (3, 4),
            'region2': (5, 6),
            'region3': (7, 8),
            'region4': (9, 10),
            'region5': (11, 12),
            'region6': (13, 14),
            'region7': (15, 16),
            'region8': (17, 18)
        }[region]

        with open(path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            header = next(reader)
            header = next(reader)
            header = next(reader)

            datas = [(float(r[column_num[0]]), float(r[column_num[1]])) for r in reader]

        if len(datas) != 8760:
            raise ValueError('読み込んだデータ数が8760個ではありません。')

        datas2 = []
        if reading_order:
            datas2.append(datas[8759])
            for i in range(8759):
                datas2.append(datas[i])
        else:
            for d in datas:
                datas2.append(d)

        date_rows = []
        dt = datetime.datetime(year=2015, month=1, day=1, hour=0, minute=0)

        for i in range(8760):
            date_rows.append(dt)
            dt = dt + datetime.timedelta(hours=1)

        self.date_rows = date_rows

        air_rows = []

        for i, d in enumerate(datas2):
            try:
                air = Air(Temperature(d[0], 'Celsius'), d[1] / 1000, 'absolute')
            except:
                if error_print:
                    print('行' + str(i) + ': エラー(over RH 100%)')
                air = Air(Temperature(d[0], 'Celsius'), 100.0, 'relative')
            air_rows.append(air)

        self.air_rows = air_rows

    def get_row_by_index(self, index):
        return (self.date_rows[index], self.air_rows[index])

    def get_air_row(self):
        return self.air_rows

