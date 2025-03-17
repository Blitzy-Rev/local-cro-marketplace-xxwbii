"""
Email utility module for the Molecular Data Management and CRO Integration Platform.

This module provides functions for sending emails, including verification emails,
password reset emails, welcome emails, and various notification emails. It implements
a fully local email sending solution without external dependencies.
"""

import os  # standard library
import smtplib  # standard library
import ssl  # standard library
from email.mime.text import MIMEText  # standard library
from email.mime.multipart import MIMEMultipart  # standard library
from typing import Dict, Any, Optional, List  # standard library

from jinja2 import Template, FileSystemLoader, Environment  # v3.1+

from ..core.config import get_settings
from ..logging_config import logger
from ..constants import PROJECT_NAME

# Define the directory where email templates are stored
TEMPLATES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../templates/email')

# Set up Jinja2 environment
jinja_env = Environment(loader=FileSystemLoader(TEMPLATES_DIR))


def get_email_settings() -> Dict[str, Any]:
    """
    Gets email configuration settings from application settings.
    
    Returns:
        Dict[str, Any]: Dictionary containing email configuration settings
    """
    settings = get_settings()
    
    # Default values if settings don't exist
    email_settings = {
        "SMTP_SERVER": getattr(settings, "SMTP_SERVER", "localhost"),
        "SMTP_PORT": getattr(settings, "SMTP_PORT", 25),
        "SMTP_USERNAME": getattr(settings, "SMTP_USERNAME", ""),
        "SMTP_PASSWORD": getattr(settings, "SMTP_PASSWORD", ""),
        "SMTP_USE_TLS": getattr(settings, "SMTP_USE_TLS", False),
        "SMTP_FROM_EMAIL": getattr(settings, "SMTP_FROM_EMAIL", "noreply@example.com"),
        "SMTP_FROM_NAME": getattr(settings, "SMTP_FROM_NAME", PROJECT_NAME),
        "SMTP_ENABLED": getattr(settings, "SMTP_ENABLED", False),
    }
    
    return email_settings


def send_email(to_email: str, subject: str, html_content: str, text_content: str = None) -> bool:
    """
    Sends an email using SMTP.
    
    Args:
        to_email: Recipient email address
        subject: Email subject
        html_content: HTML content of the email
        text_content: Plain text content of the email (optional)
        
    Returns:
        bool: True if email was sent successfully, False otherwise
    """
    email_settings = get_email_settings()
    
    # Check if email is enabled
    if not email_settings.get("SMTP_ENABLED", False):
        logger.warning("Email sending is disabled in settings")
        return False
    
    # Create message
    msg = MIMEMultipart("alternative")
    msg["Subject"] = subject
    msg["From"] = f"{email_settings['SMTP_FROM_NAME']} <{email_settings['SMTP_FROM_EMAIL']}>"
    msg["To"] = to_email
    
    # Attach text part if provided
    if text_content:
        part1 = MIMEText(text_content, "plain")
        msg.attach(part1)
    
    # Attach HTML part
    part2 = MIMEText(html_content, "html")
    msg.attach(part2)
    
    try:
        # Create secure connection with SMTP server and send email
        context = ssl.create_default_context() if email_settings["SMTP_USE_TLS"] else None
        
        with smtplib.SMTP(email_settings["SMTP_SERVER"], email_settings["SMTP_PORT"]) as server:
            if email_settings["SMTP_USE_TLS"]:
                server.starttls(context=context)
            
            # Login if credentials are provided
            if email_settings["SMTP_USERNAME"] and email_settings["SMTP_PASSWORD"]:
                server.login(email_settings["SMTP_USERNAME"], email_settings["SMTP_PASSWORD"])
            
            server.send_message(msg)
        
        logger.info(f"Email sent successfully to {to_email}", {"subject": subject})
        return True
    except Exception as e:
        logger.error(f"Failed to send email to {to_email}: {str(e)}", 
                    {"error": str(e), "subject": subject})
        return False


def render_template(template_name: str, context: Dict[str, Any]) -> str:
    """
    Renders an email template with provided context.
    
    Args:
        template_name: Name of the template file
        context: Dictionary containing template variables
        
    Returns:
        str: Rendered HTML content
    """
    # Add project name to context
    context["project_name"] = PROJECT_NAME
    
    # Get template and render with context
    template = jinja_env.get_template(template_name)
    return template.render(**context)


def send_verification_email(email: str, username: str, token: str, verification_url: str) -> bool:
    """
    Sends an email verification link to a user.
    
    Args:
        email: Recipient email address
        username: Username of the recipient
        token: Verification token
        verification_url: Base URL for verification
        
    Returns:
        bool: True if email was sent successfully, False otherwise
    """
    context = {
        "username": username,
        "verification_url": f"{verification_url}?token={token}",
        "token": token
    }
    
    html_content = render_template("verification_email.html", context)
    subject = f"Verify your email for {PROJECT_NAME}"
    
    return send_email(email, subject, html_content)


def send_password_reset_email(email: str, username: str, token: str, reset_url: str) -> bool:
    """
    Sends a password reset link to a user.
    
    Args:
        email: Recipient email address
        username: Username of the recipient
        token: Password reset token
        reset_url: Base URL for password reset
        
    Returns:
        bool: True if email was sent successfully, False otherwise
    """
    context = {
        "username": username,
        "reset_url": f"{reset_url}?token={token}",
        "token": token
    }
    
    html_content = render_template("password_reset.html", context)
    subject = f"Password Reset for {PROJECT_NAME}"
    
    return send_email(email, subject, html_content)


def send_welcome_email(email: str, username: str) -> bool:
    """
    Sends a welcome email to a newly verified user.
    
    Args:
        email: Recipient email address
        username: Username of the recipient
        
    Returns:
        bool: True if email was sent successfully, False otherwise
    """
    context = {
        "username": username
    }
    
    html_content = render_template("welcome_email.html", context)
    subject = f"Welcome to {PROJECT_NAME}"
    
    return send_email(email, subject, html_content)


def send_experiment_status_notification(
    email: str, 
    username: str, 
    experiment_name: str, 
    status: str, 
    details: str = None
) -> bool:
    """
    Sends an email notification about experiment status change.
    
    Args:
        email: Recipient email address
        username: Username of the recipient
        experiment_name: Name of the experiment
        status: New status of the experiment
        details: Additional details (optional)
        
    Returns:
        bool: True if email was sent successfully, False otherwise
    """
    context = {
        "username": username,
        "experiment_name": experiment_name,
        "status": status,
        "details": details
    }
    
    html_content = render_template("experiment_status.html", context)
    subject = f"Experiment Status Update: {experiment_name} is now {status}"
    
    return send_email(email, subject, html_content)


def send_submission_notification(
    email: str, 
    cro_name: str, 
    experiment_name: str, 
    pharma_company: str, 
    submission_details: Dict[str, Any]
) -> bool:
    """
    Sends an email notification about new submission to CRO.
    
    Args:
        email: Recipient email address (CRO)
        cro_name: Name of the CRO
        experiment_name: Name of the experiment
        pharma_company: Name of the pharmaceutical company
        submission_details: Dictionary containing submission details
        
    Returns:
        bool: True if email was sent successfully, False otherwise
    """
    context = {
        "cro_name": cro_name,
        "experiment_name": experiment_name,
        "pharma_company": pharma_company,
        "submission_details": submission_details
    }
    
    html_content = render_template("submission_notification.html", context)
    subject = f"New Experiment Submission: {experiment_name} from {pharma_company}"
    
    return send_email(email, subject, html_content)


def send_quote_notification(
    email: str, 
    username: str, 
    experiment_name: str, 
    cro_name: str, 
    price: float, 
    turnaround_days: int
) -> bool:
    """
    Sends an email notification about quote provided by CRO.
    
    Args:
        email: Recipient email address (pharma user)
        username: Username of the recipient
        experiment_name: Name of the experiment
        cro_name: Name of the CRO
        price: Quoted price
        turnaround_days: Estimated turnaround time in days
        
    Returns:
        bool: True if email was sent successfully, False otherwise
    """
    context = {
        "username": username,
        "experiment_name": experiment_name,
        "cro_name": cro_name,
        "price": price,
        "turnaround_days": turnaround_days
    }
    
    html_content = render_template("quote_notification.html", context)
    subject = f"Quote Received for {experiment_name} from {cro_name}"
    
    return send_email(email, subject, html_content)


def send_results_notification(
    email: str, 
    username: str, 
    experiment_name: str, 
    cro_name: str, 
    result_url: str = None
) -> bool:
    """
    Sends an email notification about results uploaded by CRO.
    
    Args:
        email: Recipient email address (pharma user)
        username: Username of the recipient
        experiment_name: Name of the experiment
        cro_name: Name of the CRO
        result_url: URL to view results (optional)
        
    Returns:
        bool: True if email was sent successfully, False otherwise
    """
    context = {
        "username": username,
        "experiment_name": experiment_name,
        "cro_name": cro_name,
        "result_url": result_url
    }
    
    html_content = render_template("results_notification.html", context)
    subject = f"Results Available for {experiment_name} from {cro_name}"
    
    return send_email(email, subject, html_content)


def send_system_notification(emails: List[str], subject: str, message: str) -> Dict[str, Any]:
    """
    Sends a system notification email to users.
    
    Args:
        emails: List of recipient email addresses
        subject: Email subject
        message: Notification message
        
    Returns:
        Dict[str, Any]: Dictionary with success count and total count
    """
    success_count = 0
    context = {
        "message": message
    }
    
    html_content = render_template("system_notification.html", context)
    
    for email in emails:
        if send_email(email, subject, html_content):
            success_count += 1
    
    return {
        "success_count": success_count,
        "total_count": len(emails)
    }


def is_email_enabled() -> bool:
    """
    Checks if email functionality is enabled in the application settings.
    
    Returns:
        bool: True if email is enabled, False otherwise
    """
    email_settings = get_email_settings()
    return email_settings.get("SMTP_ENABLED", False)


class EmailTemplate:
    """
    Class for managing email templates.
    """
    
    def __init__(self):
        """
        Initializes the EmailTemplate class.
        """
        self._templates = {}
        # Load default templates
    
    def load_template(self, template_name: str) -> Template:
        """
        Loads a template from file.
        
        Args:
            template_name: Name of the template file
            
        Returns:
            Template: Loaded Jinja2 template
        """
        if template_name not in self._templates:
            self._templates[template_name] = jinja_env.get_template(template_name)
        
        return self._templates[template_name]
    
    def render(self, template_name: str, context: Dict[str, Any]) -> str:
        """
        Renders a template with provided context.
        
        Args:
            template_name: Name of the template file
            context: Dictionary containing template variables
            
        Returns:
            str: Rendered template content
        """
        # Add project name to context
        context["project_name"] = PROJECT_NAME
        
        # Get template and render
        template = self.load_template(template_name)
        return template.render(**context)