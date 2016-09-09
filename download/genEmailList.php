
<html><head></head><body>

<?php
$adminPassword = "goG0nfs!m";
$word = trim($_POST['theword']);
if(strcmp($word,$adminPassword) == 0) {
	# if password was correct, login--new database.
	mysql_connect('mysql.emonet.biology.yale.edu','nfwebreguser','nepo4nEFsim') or die('Error 1.');
	@mysql_select_db('nfsimwebregister2') or die('Error 2.');

	#Get a list of all the users and information
	$getUserQuery="SELECT * FROM nfRegUser;";
	$userList=mysql_query($getUserQuery);
	$rows=mysql_numrows($userList);

	$i=0;
	//echo("<br><br>\n<table border=1 cellpadding=6>\n");
	//echo("<tr><td>ID&nbsp</td><td>NAME</td><td>EMAIL</td><td>PASSWORD</td><td>AFFILIATION</td><td>DOWNLOADS</td><td>DATE CREATED</td></tr>\n");
	while($i<$rows) {
		$id = mysql_result($userList,$i,"id");
		$name = mysql_result($userList,$i,"name");
		$email = mysql_result($userList,$i,"email");
		$psswd = mysql_result($userList,$i,"passwd");
		$aff = mysql_result($userList,$i,"affiliation");
		$dlc = mysql_result($userList,$i,"downLoadCount");
		$created = mysql_result($userList,$i,"dateRegistered");
	
		echo("$name, $email<br>");
		//echo("<tr><td>$id</td><td>$name</td><td>$email</td><td>$psswd</td><td>$aff</td><td>$dlc</td><td>$created</td></tr>\n");
		$i++;	
	};
	mysql_close();
	//echo("</table>\n<br><br><br>\n");
} else {
?>
<form name="password" method="post" action="">
<input type="text" name="theword" size="70" rows="1"></textarea>
<input type="submit" name="submit" value="+" />
<?php } ?>

</body></html>