import random
import string


class GenerateJobId(object):
    def __init__(self, lenx=5):
        """
            Generates a random job id for the job

            :param lenx: The length of the generated id
            :type lenx: int
        """

        self.lenx = lenx

    def generate(self):
        chars = string.ascii_lowercase

        id = ''.join(
            random.choice(chars) for _ in range(self.lenx)
        )

        return id
