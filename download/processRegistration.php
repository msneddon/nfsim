<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">

<html>
<head>
	<title>NFsim - register for download</title>

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
				<a href="../pages/features.html" class="topMenu">features</a> &nbsp&nbsp&nbsp
				<a href="../pages/feedback.php" class="topMenu">feedback</a> <br>
				<a href="../pages/support/support.html" class="topMenu">support</a> &nbsp&nbsp&nbsp
				<a href="../pages/models/models.html" class="topMenu">models</a> &nbsp&nbsp&nbsp
				<a href="../pages/developers.html" class="topMenu">developers</a> &nbsp&nbsp&nbsp
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
		<br>
		<div class="sectionTitleDiv">register to download NFsim</div>
		</td>
	</tr>
	
	
	<tr>
		<td colspan=4>
		<center>
		<!--  CONTENT -->

		<table width=650> <tr><td>


<?php

	//first check if anything was posted, and if so, if
	//it seems to be correct...
	$name = trim($_POST['name']);
	$email = trim($_POST['email']);
	$aff = trim($_POST['affiliation']);
	$theHow = trim($_POST['theHow']);
	$theWhy = trim($_POST['theWhy']);
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
		
	} else if(empty($aff)) {
		echo("<center><br><br>\n");
		echo("<font color=#FF0000>You did not enter your affiliation!</font><br><br>\n");
		echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
		echo("</center>\n");
		$isGood = 0;	
		
	} else if(empty($theWhy)) {
		echo("<center><br><br>\n");
		echo("<font color=#FF0000>You did not enter why you want to use NFsim!</font><br><br>\n");
		echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
		echo("</center>\n");
		$isGood = 0;	
		
	}
	
	
	//Check if the email address is valid
	if($isGood==1) {
		$containsSpace = stripos($email," ");
		$atSymbolPos = stripos($email,"@");
		if(!empty($containsSpace) or empty($atSymbolPos)) {
			echo("<center><br><br>\n");
			echo("<font color=#FF0000>The email address you entered ($email) appears to be invalid!</font><br><br>\n");
			echo("We need a valid email to provide you with the download information!<br><br>\n");
			echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
			echo("</center>\n");
			$isGood = 0;
		}
	}
	
	
	

	
	//Check the captcha

	/*if($isGood==1) {
		require_once('recaptchalib.php');
		$privatekey = "6LfOzAkAAAAAAFXOEbsBA74XpVoOx5m-amESq11I";
		//$privatekey = "6Let67kSAAAAAJ5rsEKfLpkczoaW-1OrrhHBOe-D";
		
 	    $resp = null;
 	    $error = null;
 	    
 	    
 	    $resp = recaptcha_check_answer ($privatekey,
                                        $_SERVER["REMOTE_ADDR"],
                                        $_POST["recaptcha_challenge_field"],
                                        $_POST["recaptcha_response_field"]);
                                        
		if (!$resp->is_valid) {
			echo("<center><br><br>\n");
			echo("<font color=#FF0000>Sorry! You did not correctly solve the CAPTCHA puzzle.</font><br><br>\n");
			echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
			echo("The error given by reCAPTCHA was: <font color=#FF0000>$resp->error</font><br>\n");
			echo("</center>\n");
			$isGood=0;
		}
		//else {
		//	echo("good\n");
		//}
		

	}*/

	//srinivas-disabled catpchpa everywhere for testing
	


	
	// log on to the database to begin our operations...
	//mysql_connect('localhost','nfWebRegUser','nepo4nEFsim') or die('internal database error 1.');
	//@mysql_select_db('NFsimWebRegister') or die('internal database error 2.');

	// new code to log into new database...
	mysql_connect('mysql.emonet.biology.yale.edu','nfwebreguser','nepo4nEFsim') or die('internal database error 1.');
	@mysql_select_db('nfsimwebregister2') or die('internal database error 2.');


	
	//Got here, so everything appears to be in order
	if($isGood==1) {
		
		// first check if the name exists or not
		$getUserQuery="SELECT id,email FROM nfRegUser;";
		$userList=mysql_query($getUserQuery);
		$rows=mysql_numrows($userList);
		$i=0;
		while($i<$rows) {
			$cEmail = mysql_result($userList,$i,"email");
			if(strcmp($cEmail,$email) == 0) {
				echo("<center><br><br>\n");
				echo("<font color=#FF0000>The email you provided has already been used to register for NFsim!.</font><br><br>\n");
				echo("If you forgot your password, or think this is an error, email <a href=\"mailto:mwsneddon@lbl.gov\">mwsneddon@lbl.gov</a>.<br><br>\n");
				echo("Please go back by pressing the back button on your browser and try again.<br><br>\n");
				echo("</center>\n");
				$isGood=0;
				break;
			}
			$i++;
		}
	}

	if($isGood==1) {
		// insert the entry into the database
		#$randNum = rand(100,999);
		$randNum = substr(md5(rand().rand()), 0, 6);
		$autoGenPassword = "goNF_".$randNum;
		#echo("password = $autoGenPassword <br>");
		
		//make sure we don't have any quotes! This will mess up our sql insert!
		$aff = str_replace("'"," ",$aff);
		$aff = str_replace("\""," ",$aff);
		
		$theHow = str_replace("'"," ",$theHow);
		$theHow = str_replace("\""," ",$theHow);
		
		$theWhy = str_replace("'"," ",$theWhy);
		$theWhy = str_replace("\""," ",$theWhy);
		
		
		$insertStatement="INSERT INTO nfRegUser (name, email, passwd, affiliation, theHow, theWhy, downLoadCount)";
		$insertStatement=$insertStatement."VALUES ('$name','$email','$autoGenPassword','$aff','$theHow','$theWhy',0)";
		//echo("Insert: $insertStatement <br>");
		
		if (!mysql_query($insertStatement))  {
  			die('internal database error 3.  This error may be caused if you included a single or double quote or apostrophe in one of your forms.');
  		}
	}
		
	// clse the door
	mysql_close();
	
	
	//If things are good, send out an email
	if($isGood==1) {
		
		// send one to us admin folk
		$adminMssg = "\nUser Information\n--------------------------\n";
		$adminMssg = $adminMssg."Name: ".$name."\n";
		$adminMssg = $adminMssg."Email: ".$email."\n";
		$adminMssg = $adminMssg."Password: ".$autoGenPassword."\n";
		$adminMssg = $adminMssg."Affiliation: ".$aff."\n";
		$adminMssg = $adminMssg."The How: ".$theHow."\n";
		$adminMssg = $adminMssg."The Why: ".$theWhy."\n";
		
		$adminSubj = "New NFsim User Registered!";
		
		require_once('emailAlerts.php');
		sendEmailMessage("mwsneddon@lbl.gov", $adminSubj, $adminMssg);
		sendEmailMessage("nfsim.super.user@gmail.com", $adminSubj, $adminMssg);
		sendEmailMessage("thierry.emonet@yale.edu", $adminSubj, $adminMssg);
		sendEmailMessage("faeder@pitt.edu", $adminSubj, $adminMssg);  
		sendEmailMessage("emonetwebsite@srinivas.gs", $adminSubj, $adminMssg);
		
		
		// Send one to the new user
		$userMssg = "Hello $name,\n\n";
		$userMssg = $userMssg."Thank you for registering! To download NFsim now, go to the\n";
		$userMssg = $userMssg."download page: http://emonet.biology.yale.edu/nfsim/download/download.php\n";
		$userMssg = $userMssg."\nWhen you get there, enter the email address you registered\n";
		$userMssg = $userMssg."with and your password. Your password is: $autoGenPassword\n\n";
		$userMssg = $userMssg."Happy simulating!\n";
		$userSubj = "Welcome to NFsim!";
		
		sendEmailMessage($email, $userSubj, $userMssg);
	}
		
		
		
	if($isGood==1) {	
		// give the user a nice message
		echo("<br><br>");
		echo("Thank you for your interest in NFsim!  Your registration was successful!<br><br>\n");
		echo("You will recieve an email in a few minutes with instructions and a password for downloading\n");
		echo("NFsim. <br><br> \n");
		echo("If you do not recieve the email, please contact <a href=\"mailto:mwsneddon@lbl.gov\">mwsneddon@lbl.gov</a>\n");
		echo("as your email address may have been incorrectly processed.<br><br>");
		
		echo("<br><br>");
		echo("Here is the information you submitted:<br><br>");
		echo("&nbsp&nbsp&nbsp&nbsp&nbsp <b>$name</b><br>");
		echo("&nbsp&nbsp&nbsp&nbsp&nbsp <i>$aff</i><br><br>");
		echo("&nbsp&nbsp&nbsp&nbsp&nbsp $email<br>");
		
	}

	
?>








</td></tr></table>
<!-- END CONTENT -->
		
	</td>
	</tr>
</table>
</center>

</body>
</html>