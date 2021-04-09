# Here lies the Celery app import
from ..celery_app import app

# Here we import the Celery Task class
from celery import Task

# Here we import Hunana ONLY FOR TYPING
from hunana.datatypes import Position as HunanaPosition
from hunana.datatypes import Variant as HunanaVariant

# Here lies the model imports
from ..models import HunanaPosition, Job, Result

# Here lies other imports
from .logging import Logging

# Here lies constants imports
from .constants import LogContexts

# Here lies Python default library imports
from itertools import product
from collections import namedtuple


class ATAT(Task):
    name = "ATAT"
    queue = "ATAT"

    switch = namedtuple('MotifSwitch', ['position', 'sequence', 'fromx', 'to'])

    def run(self, hunana_results: list) -> list:
        """
            Stores the Hunana results for Source and Reservoir under the appropriate job id.

            :param hunana_results: A list containing two dicts pertaining to the tagged Hunana results for Source a
            Reservoir sequences.
            :type hunana_results: list

            :returns: A list of NamedTuple pertaining to each switch observed from Source --> Reservoir
        """

        Logging.make_log_entry(LogContexts.INFO, 'Starting transmissibility analysis.')

        positions = zip(*hunana_results)
        switches = list()

        for source_position, reservoir_position in positions:  # type: dict
            for source_variant, reservoir_variant in product(source_position.get('variants'),
                                                             reservoir_position.get('variants')):  # type: dict
                if source_variant.get('sequence') == reservoir_variant.get('sequence'):
                    if source_variant.get('motif_short') != reservoir_variant.get('motif_short'):
                        switches.append(
                            {
                                'position': source_variant.get('position'),
                                'sequence': source_variant.get('sequence'),
                                'from': source_variant.get('motif_long').upper(),
                                'to': reservoir_variant.get('motif_long').upper()
                            }
                        )

        Logging.make_log_entry(LogContexts.INFO, 'Transmissibility analysis completed.')

        return switches


@app.task(name="Noop")
def noop(hunana_results: list) -> list:
    """
        This performs no actions. Just returns the Hunana results so we can inject those results into the group
    """

    return hunana_results


# Here we register the task under our Celery app  because when ug Class based tasks auto-register does not happen
app.register_task(ATAT())
