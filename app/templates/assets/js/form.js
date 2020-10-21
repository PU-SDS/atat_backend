const JOB_FORM_ID = "submit-job";
const SUBMIT_JOB_BUTTON_ID = "submit-job-btn"

$(function() {
    $("#" + JOB_FORM_ID).submit(function(e) {
        e.preventDefault();

        $("#" + SUBMIT_JOB_BUTTON_ID).prop("disabled", true);

        var endpoint = e.currentTarget.action;
        var form_data = new FormData(e.currentTarget);

        $.ajax({
            url: endpoint,
            type: "POST",
            enctype: "multipart/form-data",
            data: form_data,
            processData: false,
            contentType: false,
            cache: false,
            timeout: 600000
        }).done(function(result) {
            $("#" + SUBMIT_JOB_BUTTON_ID).prop("disabled", false);
        })
    })
})