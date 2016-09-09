<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">

<html>
<head>
	<title>NFsim - feedback</title>

	<META NAME="description" CONTENT=".">
	<META NAME="keywords" CONTENT="">
	<link rel="stylesheet" type="text/css" href="../css/nf_style.css" />

</head>

<body>


	<center>
	<table border=0 cellspacing=0 cellpadding=0>
	<tr>
		<!-- <td width=30%></td> -->
		<td><img src="../images/topLogo.jpg" align=right></td>
		<td bgcolor=#2E3192 valign="center">
			<center>
			<div class="topMenuDiv">
				<a href="../" class="topMenu">overview</a>&nbsp&nbsp&nbsp
				<a href="../download" class="topMenu">download</a>&nbsp&nbsp&nbsp
				<a href="features.html" class="topMenu">features</a> &nbsp&nbsp&nbsp
				<a href="feedback.html" class="topMenu">feedback</a> <br>
				<a href="support/support.html" class="topMenu">support</a> &nbsp&nbsp&nbsp
				<a href="models/models.html" class="topMenu">models</a> &nbsp&nbsp&nbsp
				<a href="developers.html" class="topMenu">developers</a> &nbsp&nbsp&nbsp
				<a href="http://emonet.biology.yale.edu" class="topMenu">emonet lab</a>
			</div>
			</center>
		</td>
		<td bgcolor=#2E3192 width=15%></td>
		<td><img src="../images/topRight.jpg" align=left></td>
		<!-- <td bgcolor=#2E3192 width=50%></td> -->
	</tr>
	
	<tr>
		<td colspan=4>
		
		<div class="sectionTitleDiv">
			feedback</div>
		</td>
	</tr>
	
	
	<tr>
		<td colspan=4>
		<center>
		<!--  CONTENT -->
		
		
		<table border=0 width=650><tr><td>
		
		
		
		<?php

	//first check if anything was posted, and if so, if
	//it seems to be correct...
	$name = trim($_POST['name']);
	$email = trim($_POST['email']);
	$feedback = trim($_POST['feedback']);
	
	
	if(empty($name)) {
		
?>

		<p>Please submit your questions, comments, or other requests here.  We will try to
		respond as quickly as possible, and may post your question and our
		response to the <a href="support/support.html">support</a> page to help other users.
		</p>
		
		
		
		</td></tr>
		</table>
		
		<form name="registerForm" method="post" action="feedback.php">

			<center>
			
				<table border=0 cellspacing=1 cellpadding=4>
					<tr><td>
						<div class="regTextDiv">Name<font color=#FF0000>*</font></div></td>
						<td><input type="text" name="name" size="60" rows="1"></textarea> </td></tr>
					
					<tr><td>
						<div class="regTextDiv">Return E-mail<font color=#FF0000>*</font></div></td>
						<td><input type="text" name="email" size="60" rows="1"></textarea> </td></tr>
						
					<tr><td>
						<div class="regTextDiv">Question / Comment<font color=#FF0000>*</font></div></td>
						<td><textarea type="text" name="feedback" cols="50" rows="15"></textarea></td></tr>
					
					<tr><td></td><td>
						<div class="regTextRequiredFieldDiv">*Required fields</td></tr>
				</table>
	
				<br><br>
					<i>
					<?php
					//require_once('../download/recaptchalib.php');
					//$publickey = "6LfOzAkAAAAAAOv25glya8vNLI3ifEZY27rFcQLo";
					//echo recaptcha_get_html($publickey);
					?>
					<br/>
					<input type="submit" name="submit" value="Send Feedback" />
			</center>
		</form>
		
		
<?php
	} else {
	
		$isGood = 1;
		
		if(empty($name)) {
			echo("<center><br><br>\n");
			echo("<font color=#FF0000>You did not enter your name!</font><br><br>\n");
			echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
			echo("</center>\n");
			$isGood = 0;
			
		} else if(empty($email)) {
			echo("<center><br><br>\n");
			echo("<font color=#FF0000>You did not enter your email address!</font><br><br>\n");
			echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
			echo("</center>\n");
			$isGood = 0;
			
		} else if(empty($feedback)) {
			echo("<center><br><br>\n");
			echo("<font color=#FF0000>You did not enter any feedback message!</font><br><br>\n");
			echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
			echo("</center>\n");
			$isGood = 0;	
			
		}
		
		//Check the captcha
		// if($isGood==1) {
		//	require_once('../download/recaptchalib.php');
		//	$privatekey = "6LfOzAkAAAAAAFXOEbsBA74XpVoOx5m-amESq11I";
		//	$resp = recaptcha_check_answer ($privatekey,
	   //                             $_SERVER["REMOTE_ADDR"],
	   //                             $_POST["recaptcha_challenge_field"],
	  //                              $_POST["recaptcha_response_field"]);
	
		//	if (!$resp->is_valid) {
		//		echo("<center><br><br>\n");
		//		echo("<font color=#FF0000>Sorry! You did not correctly solve the CAPTCHA puzzle.</font><br><br>\n");
		//		echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
		//		echo("The error given by reCAPTCHA was: <font color=#FF0000>$resp->error</font><br>\n");
		//		echo("</center>\n");
		//		$isGood=0;
		//	}
		//}
		
		
		//Check if the email address is valid
		if($isGood==1) {
			$containsSpace = stripos($email," ");
			$atSymbolPos = stripos($email,"@");
			if(!empty($containsSpace) or empty($atSymbolPos)) {
				echo("<center><br><br>\n");
				echo("<font color=#FF0000>The email address you entered ($email) appears to be invalid!</font><br><br>\n");
				echo("We need a valid email in order for us to respond to your request!<br><br>\n");
				echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
				echo("</center>\n");
				$isGood = 0;
			}
		}
		
		//If things are good, send out the feedback email
		if($isGood==1) {
			
			// send one to us admin folk
			$adminMssg = "\n";
			$adminMssg = $adminMssg."Sent from: '".$name."'\n";
			$adminMssg = $adminMssg."Return address: '".$email."'\n\n";
			$adminMssg = $adminMssg."Feedback:\n--------------------------\n";
			$adminMssg = $adminMssg.$feedback."\n";
			
			$adminSubj = "New NFsim Feedback Message Sent!";
			
			require_once('../download/emailAlerts.php');
			sendEmailMessage("mwsneddon@lbl.gov", $adminSubj, $adminMssg);
		}
			
			
		if($isGood==1) {	
			// give the user a nice message
			echo("<br><br><center>");
			echo("Thank you for your feedback!  Your message was sent successfully!</center><br><br><br><br>\n");
			
		}
		// else {
			// give the user a nice message
		//	echo("<br><br><center>");
		//	echo("There was an error processing your message.<br>Please press the back button on your browser and try again.</center><br><br>\n");
		//}

	}

	
?>
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		</td></tr>
		</table>
		

		
		<br><br><br><br><br><br><br><br>
		
		
	</td>
	</tr>

	
	
	

</table>
</center>

</body>
</html>
