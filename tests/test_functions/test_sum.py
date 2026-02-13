import main

def test_calc_raw_tax_set_1(test1):
    res = main.calc_raw_tax_set(test1["file"]["header_line"], test1["file"]["lines"])
    assert res[2] == False