OK = "ok ✅"
NOK = "❌"
EXCL = "❗"
NOTHING_TO_DO = "💤"
HURRAH = "🥳"

class ANSI(object):
    """ Class defining some ANSI control sequences, for example for
    changing text color """
    CSI = '\x1b[' # ANSI Control Sequence Introducer. Not the TV show
    OSC = '\x1b]8'
    SGR = 'm' # Set Graphics Rendition code
    NORMAL = CSI + '0' + SGR
    BLACK = CSI + '30' + SGR
    GREY = CSI + '30;1' + SGR
    RED = CSI + '31' + SGR
    BRIGHTRED = CSI + '31;1' + SGR
    GREEN = CSI + '32' + SGR
    BRIGHTGREEN = CSI + '32;1' + SGR
    YELLOW = CSI + '33' + SGR
    BRIGHTYELLOW = CSI + '33;1' + SGR
    BLUE = CSI + '34' + SGR
    BRIGHTBLUE = CSI + '34;1' + SGR
    MAGENTA = CSI + '35' + SGR
    BRIGHTMAGENTA = CSI + '35;1' + SGR
    CYAN = CSI + '36' + SGR
    BRIGHTCYAN = CSI + '36;1' + SGR
    # There is another grey with code "37" - white without intensity
    # Not sure if it is any different from "30;1" aka "bright black"
    WHITE = CSI + '37;1' + SGR
    def hyperlink(target, text=None, width=None):
        if text is None:
            text = target
        text = str(text)
        pad = ''
        if width is not None:
            text = text[:width]
            pad = ' ' * (width - len(text))
        return "{0};;{1}\x07{2}{0};;\x07".format(ANSI.OSC, target, text) + pad



def print_command_line(what, *args):
    line = " ".join([what] + list(args))
    print(f"{ANSI.YELLOW}{what} command line: {line}{ANSI.NORMAL}")

def major_message(msg, *args, **kwargs):
    print(ANSI.BRIGHTBLUE + str(msg), *args, ANSI.NORMAL, **kwargs)

def warning_message(*args, **kwargs):
    print(ANSI.RED + "Warning:", *args, ANSI.NORMAL, **kwargs)

def error_message(*args, **kwargs):
    print(ANSI.BRIGHTRED + "Error:", *args, ANSI.NORMAL, **kwargs)

