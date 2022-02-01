from Moore import options
from Hlt2Conf.lines.b_to_kpi import all_lines

def make_lines():
    return [builder() for builder in all_lines.values()]

options.lines_maker = make_lines
