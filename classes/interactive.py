def answer_yes(input_key):
    return input_key.upper() == "Y" or input_key.upper() == "YES"


def answer_no(input_key):
    return input_key.upper() == "N" or input_key.upper() == "NO"


def answer_empty(input_key):
    return len(input_key) == 0


def ask_to_continue_input():
    result = False
    while True:
        input_key = input("Your answer is wrong. Do you want to try again? (Y[es]/n[o])> ")
        if answer_yes(input_key) or answer_empty(input_key):
            result = True
            break
        elif answer_no(input_key):
            print("See you later")
            exit()
    return result