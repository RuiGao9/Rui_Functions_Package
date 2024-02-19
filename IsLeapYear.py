def IsLeapYear(Year):
    """ If year is a leap year return True;
        else return False.
        Copy and past from:
        https://stackoverflow.com/questions/620305/convert-year-month-day-to-day-of-year-in-python
    """
    if Year % 100 == 0:
        return Year % 400 == 0
    return Year % 4 == 0