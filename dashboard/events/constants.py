from app.warehousing import JobStatuses


class PrecheckStatusStyles(object):
    STATUS_BASECLASS = 'alert alert-%color% d-flex justify-content-center'

    def set(self, status: str) -> str:
        if status == JobStatuses.JOB_NONEXIST:
            return self.STATUS_BASECLASS.replace('%color%', 'warning')
        elif status == JobStatuses.JOB_PROCESSING:
            return self.STATUS_BASECLASS.replace('%color%', 'info')
        elif status == JobStatuses.JOB_FAILED:
            return self.STATUS_BASECLASS.replace('%color%', 'danger')
        elif status == JobStatuses.JOB_COMPLETED:
            return self.STATUS_BASECLASS.replace('%color%', 'success')
