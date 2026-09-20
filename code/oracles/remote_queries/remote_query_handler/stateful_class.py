import json

from logging import Logger

class StatefulClassException(Exception):
    pass

class StatefulClass:
    def __init__(self, work_dir: str, logging: Logger,
                 attrs: list, state_fname: str, *args):
        
        self.work_dir   = work_dir
        self.logging    = logging
        self.attrs      = attrs
        self.state_file = f"{self.work_dir}/{state_fname}"
        
        if len(args) == 0:
            with open(self.state_file, "r") as fp:
                state_dict = json.load(fp)
        elif len(args) == len(attrs):
            state_dict = dict()
            # Assume argument order
            for i, attr in enumerate(self.attrs):
                state_dict[attr] = args[i]
        else:
            raise StatefulClassException(f"Invalid arguments for stateful class, expected: {', '.join(attrs)} but only got {len(args)} arguments.")

        for attr in self.attrs:
            if attr in state_dict:
                setattr(self, attr, state_dict[attr])
            else:
                self.logging.warning(f"Unknown attribute '{attr}' from {self.state_file} not in attribute list {self.attrs}.")
        self.save()
                
    def save(self):
        state_dict = dict()
        for attr in self.attrs:
            state_dict[attr] = getattr(self, attr)

        with open(self.state_file, "w") as fp:
            json.dump(state_dict, fp)
