$(document).ready(function () {

    $("#btnSubmit").click(function (event) {

        //stop submit the form, we will post it manually.
        event.preventDefault();

        // Get form
        var form = $('#atatJobSubmissionForm')[0];

		// Create an FormData object
        var data = new FormData(form);

		// If you want to add an extra field for the FormData
        //data.append("CustomField", "This is some extra data, testing");

		// disabled the submit button
        $("#btnSubmit").prop("disabled", true);

        $.ajax({
            type: "POST",
            enctype: 'multipart/form-data',
            url: "/submit",
            data: data,
            processData: false,
            contentType: false,
            cache: false,
            timeout: 600000,
	}).done(function(data){
		if(data.error){
			$("#errorAlert").text(data.error).show();
                        $("#successAlert").hide();
		} else {
                	$("#successAlert").text(data.message).show();
                	$("#errorAlert").hide();
            	}
		$("#btnSubmit").prop("disabled", false);
        });

    });

});
