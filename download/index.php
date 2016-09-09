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



<center>
(Already registered? Then click <a href="download.php">here</a>.)
</center>

<br><br>

	<center>
	<table border=0 width=650> <tr><td>
	    <p>In order to download NFsim, please complete the quick registration process below.  Registration
	    allows us to track usage and notify you when important updates are released.</p>
	    
	    <p>After registering, you will recieve an email with your password and instructions for download.
	    If you do not recieve an email or have other problems with downloading NFsim, please
	    contact <a href="mailto:mwsneddon@lbl.gov">mwsneddon@lbl.gov</a>.</p>
	    
	    <p>And if you decide to use NFsim for any research application or in a research publication, 
	    please let us know!</p>

	</td></tr></table>
	<br><br>
	</center>


<form name="registerForm" method="post" action="processRegistration.php">

	<center>
	
		<table border=0 cellspacing=1 cellpadding=4>
			<tr><td>
				<div class="regTextDiv">Name<font color=#FF0000>*</font></div></td>
				<td><input type="text" name="name" size="60" rows="1"></textarea> </td></tr>
			
			<tr><td>
				<div class="regTextDiv">E-mail<font color=#FF0000>*</font></div></td>
				<td><input type="text" name="email" size="60" rows="1"></textarea> </td></tr>
				
			<tr><td>
				<div class="regTextDiv">Affiliation<font color=#FF0000>*</font></div></td>
				<td><input type="text" name="affiliation" size="60" rows="1"></textarea></td></tr>
			
			<tr><td>
				<div class="regTextDiv">How did you hear about NFsim?</div></td>
				<td><input type="text" name="theHow" size="60" rows="1"></textarea> </td></tr>
			
			<tr><td>
				<div class="regTextDiv">Why do you want to use NFsim?<font color=#FF0000>*</font></div></td>
				<td><input type="text" name="theWhy" size="60" rows="1"></textarea> </td></tr>
		
			<tr><td></td><td>
				<div class="regTextRequiredFieldDiv">*Required fields</td></tr>
		</table>
	
		<br><br>
			<i>
			<?php
			//require_once('recaptchalib.php');
			//$publickey = "6LfOzAkAAAAAAOv25glya8vNLI3ifEZY27rFcQLo";
			//$publickey = "6Let67kSAAAAAPKN_PrY90saobVT1C6icGUIIktp";
			//echo recaptcha_get_html($publickey);
			?>>
			<br/>
			<input type="submit" name="submit" value="Register" />
	</center>
</form>

<br><br><br><br><br>


<!-- END CONTENT -->
		
	</td>
	</tr>
</table>
</center>

</body>
</html>